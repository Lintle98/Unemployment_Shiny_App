# ==================================================================
# Unemployment and Deprivation in Education and Basic Services in KwaZulu-Natal:
# Associations and Spatial Dependence (Census 2011)
# STAT7006A Spatial Statistics - group R Shiny assignment
#
# To run:    open this file in RStudio and click "Run App"
# To publish: open deploy.R and click "Source"
# Data:      data/municipalities.rds (made by 01_prepare_data.R)
# ==================================================================

library(shiny)
library(sf)
library(dplyr)
library(spdep)
library(spatialreg)
library(leaflet)
library(plotly)
library(DT)
library(bslib)      # load last: spdep also has a card() function

# ===================== DATA (loaded once) =====================
muni <- readRDS("data/municipalities.rds")
if (is.na(sf::st_crs(muni))) stop("Municipal boundary CRS is missing. Set the correct CRS in data preparation.")
muni <- sf::st_transform(muni, 4326)
muni_df <- st_drop_geometry(muni)

vars <- c("Unemployment rate (%)"            = "unemployment",
          "Youth unemployment, 15-34 (%)"    = "youth_unemp",
          "Matric, aged 20+ (%)"             = "matric",
          "No schooling, aged 20+ (%)"       = "no_school",
          "Piped water inside dwelling (%)"  = "piped_water",
          "Electricity for lighting (%)"     = "electricity")
predictors <- vars[c(3, 4, 5, 6)]
lab <- function(v) names(vars)[vars == v]

# ---- Data validation (reported in the About tab; stops on fatal problems) ----
data_checks <- data.frame(
  Check = c("Municipalities (rows)", "Unique municipality codes", "Missing values in indicators",
            "Indicators within 0-100%", "Valid geometries", "Attribute rows match polygons"),
  Result = c(
    nrow(muni),
    ifelse(anyDuplicated(muni$muni_code) == 0, "Yes", "No - duplicates found"),
    sum(is.na(muni_df[, vars])),
    ifelse(all(as.matrix(muni_df[, vars]) >= 0 & as.matrix(muni_df[, vars]) <= 100, na.rm = TRUE), "Yes", "No"),
    ifelse(all(sf::st_is_valid(muni)), "Yes", "No"),
    ifelse(nrow(muni_df) == length(sf::st_geometry(muni)), "Yes", "No")))
if (anyDuplicated(muni$muni_code) > 0) stop("Duplicate municipality codes in data/municipalities.rds.")
if (sum(is.na(muni_df[, vars])) > 0) stop("Missing indicator values in data/municipalities.rds.")

# Representative points: computed in a projected CRS (UTM 36S), then returned to lon/lat
# so that knearneigh(..., longlat = TRUE) receives longitude/latitude.
centroids <- sf::st_coordinates(sf::st_transform(
  sf::st_point_on_surface(sf::st_transform(sf::st_geometry(muni), 32736)), 4326))

# Spatial weights (Chapter 3): queen, rook or k-nearest neighbours, row-standardised
make_nb <- function(type, k = 4) {
  switch(type,
         queen = poly2nb(muni, queen = TRUE),
         rook  = poly2nb(muni, queen = FALSE),
         knn   = knn2nb(knearneigh(centroids, k = k, longlat = TRUE)))
}
make_lw <- function(type, k = 4) nb2listw(make_nb(type, k), style = "W", zero.policy = TRUE)

# Getis-Ord Gi* classes (confidence levels, as in standard hot spot analysis)
gi_cols <- c("Hot spot (99%)" = "#b2182b", "Hot spot (95%)" = "#ef8a62", "Hot spot (90%)" = "#fddbc7",
             "Not significant" = "#f2f2f2",
             "Cold spot (90%)" = "#d1e5f0", "Cold spot (95%)" = "#67a9cf", "Cold spot (99%)" = "#2166ac")
w_names <- c(queen = "queen contiguity", rook = "rook contiguity", knn = "k-nearest neighbours")

fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))

# REK four-step hypothesis answer
rek_box <- function(h0, h1, stat_txt, p, alpha, reject_txt, keep_txt) {
  decision <- if (p < alpha) {
    sprintf("p = %s < α = %s, so we <b>reject H₀</b>.", fmt_p(p), alpha)
  } else {
    sprintf("p = %s ≥ α = %s, so we <b>fail to reject H₀</b>.", fmt_p(p), alpha)
  }
  conclusion <- if (p < alpha) reject_txt else keep_txt
  HTML(paste0(
    "<ol class='rek'>",
    "<li><b>Hypotheses.</b> H₀: ", h0, "<br>H₁: ", h1, "</li>",
    "<li><b>Test statistic and p-value.</b> ", stat_txt, "</li>",
    "<li><b>Decision.</b> ", decision, "</li>",
    "<li><b>Conclusion.</b> ", conclusion, "</li></ol>"))
}

# Variance inflation factors without extra packages: VIF_j = 1 / (1 - R²_j)
vif_table <- function(preds, d) {
  if (length(preds) < 2) return(data.frame(Variable = lab(preds), VIF = 1))
  data.frame(Variable = sapply(preds, lab),
             VIF = sapply(preds, function(v) {
               r2 <- summary(lm(reformulate(setdiff(preds, v), response = v), data = d))$r.squared
               1 / (1 - r2)
             }), row.names = NULL)
}

popup_html <- paste0(
  "<b style='font-size:14px'>", muni$muni_name, "</b> (", muni$muni_code, ")<br>",
  "<small>District: ", muni$district, " · Census 2011</small><hr style='margin:4px 0'>",
  "Unemployment: <b>", muni$unemployment, "%</b><br>",
  "Youth unemployment: ", muni$youth_unemp, "%<br>",
  "Matric: ", muni$matric, "%<br>",
  "No schooling: ", muni$no_school, "%<br>",
  "Piped water: ", muni$piped_water, "%<br>",
  "Electricity: ", muni$electricity, "%")

lisa_cols <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6",
               "High-Low" = "#fdae61", "Low-High" = "#b39ddb",
               "Not significant" = "#cfe3f7")   # light blue background for non-significant areas

base_map <- function() {
  leaflet(muni, options = leafletOptions(minZoom = 6)) %>%
    addProviderTiles(providers$OpenStreetMap, group = "Streets") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addLayersControl(baseGroups = c("Streets", "Satellite"),
                     options = layersControlOptions(collapsed = TRUE))
}

# Highlight is rendered as a separate outline, not a new data class.
add_selection_outline <- function(map, selected_code) {
  if (is.null(selected_code) || !nzchar(selected_code)) return(map)
  feature <- muni[as.character(muni$muni_code) == selected_code, ]
  if (nrow(feature) == 0) return(map)
  map %>% addPolygons(data = feature, fill = FALSE, color = "#132A54",
                      weight = 5, opacity = 1, options = pathOptions(clickable = FALSE))
}

choropleth <- function(v, title, selected_code = "") {
  pal <- colorBin("YlOrRd", muni[[v]], bins = 5, pretty = TRUE)
  map <- base_map() %>%
    addPolygons(layerId = as.character(muni$muni_code), fillColor = pal(muni[[v]]), fillOpacity = 0.8,
                color = "white", weight = 0.7,
                label = paste0(muni$muni_name, ": ", muni[[v]], "% (Census 2011)"),
                popup = popup_html,
                highlightOptions = highlightOptions(weight = 3, color = "#333",
                                                    bringToFront = TRUE)) %>%
    addLegend("bottomright", pal = pal, values = muni[[v]], title = paste0(title, "<br>Census 2011"),
              labFormat = labelFormat(suffix = "%"))
  add_selection_outline(map, selected_code)
}

# Static PNG maps for download (base R graphics, no extra packages)
png_map <- function(file, fill, title, legend_lab, legend_col, caption) {
  png(file, width = 2000, height = 2200, res = 250)
  op <- par(mar = c(3, 0.5, 3, 0.5)); on.exit({ par(op); dev.off() })
  plot(sf::st_geometry(muni), col = fill, border = "white", lwd = 0.6, main = title)
  legend("bottomright", legend = legend_lab, fill = legend_col, bty = "n", cex = 0.85)
  mtext(caption, side = 1, line = 1.5, cex = 0.7)
}


# Small interactive tables: sortable columns, hover highlighting, no paging for short tables
dt_small <- function(df, digits = 2) {
  df <- as.data.frame(df)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits)
  datatable(df, rownames = FALSE, class = "compact hover stripe",
            options = list(dom = if (nrow(df) > 12) "tip" else "t", pageLength = 12, ordering = TRUE))
}

# Plotly figures: high-resolution PNG from the camera button
pconf <- function(p, name) config(p, toImageButtonOptions = list(format = "png", filename = name, scale = 3))

# ===================== USER INTERFACE =====================
# Colours match the group presentation (navy, blue, teal, gold)
navy <- "#12233F"; blue <- "#3B6CF6"; purple <- "#7B5CFF"; teal <- "#16B8A6"; gold <- "#E0A526"
theme <- bs_theme(version = 5, bg = "#F7F8FA", fg = "#2B3442",
                  primary = navy, secondary = blue, info = teal, warning = gold,
                  success = teal, danger = "#C8553D",
                  base_font = font_google("Inter", local = FALSE),
                  heading_font = font_google("Inter", local = FALSE))

obj_card <- function(num, title, text, col) {
  card(class = "obj-card",
       div(class = "obj-num", style = paste0("background:", col, ";"), num),
       h5(title, class = "obj-title"), p(text, class = "mb-0"))
}
section_note <- function(...) p(class = "note", ...)

ui <- page_navbar(
  title = "KZN Unemployment & Deprivation Explorer",
  theme = theme,
  fillable = FALSE,
  bg = navy,
  header = tags$style(HTML("
    .rek li { margin-bottom: 6px; }
    .note { color: #555; font-size: 0.92em; }
    .value-box-value { font-size: 1.6rem !important; }
    .hero { background: #12233F; color: #fff; border-radius: 14px; padding: 26px 30px; margin-bottom: 18px; }
    .hero h2 { color: #fff; font-weight: 700; margin-bottom: 14px; font-size: calc(1.1rem + 0.9vw); }
    .hero .lbl { color: #E0A526; font-weight: 700; letter-spacing: 2px; font-size: 0.8rem; }
    .hero p { font-size: 1.05rem; }
    .obj-card { background: #E9EEF7; border: none; }
    .obj-num { width: 46px; height: 46px; border-radius: 50%; color: #fff; font-weight: 700;
               font-size: 1.4rem; display: flex; align-items: center; justify-content: center; margin-bottom: 8px; }
    .obj-title { color: #0F4C45; font-weight: 700; }
    .key { background: #12233F; color: #fff; border-radius: 10px; padding: 14px 18px; margin-top: 10px; }
    .key b { color: #E0A526; }
    .card-header { font-weight: 600; color: #0F4C45; }
    .glance td { padding: 2px 6px; font-size: 0.88em; }
    .scenario { background: #FFF4D6; border-radius: 8px; padding: 10px 14px; margin-bottom: 12px; }
  ")),
  sidebar = sidebar(
    width = 300,
    selectInput("var", "Variable", choices = vars),
    selectizeInput("selected_muni", "Highlight municipality (or click a map)",
                   choices = c("None" = "", stats::setNames(as.character(muni$muni_code), muni$muni_name)),
                   selected = ""),
    uiOutput("glance"),
    radioButtons("wtype", "Spatial weights",
                 c("Queen contiguity" = "queen", "Rook contiguity" = "rook",
                   "k-nearest neighbours" = "knn")),
    conditionalPanel("input.wtype == 'knn'",
                     sliderInput("k", "k (neighbours)", 2, 8, 4, step = 1)),
    radioButtons("alpha", "Significance level (α)", c("0.05" = 0.05, "0.01" = 0.01),
                 inline = TRUE),
    checkboxInput("lisa_fdr", "LISA: control false discovery rate (Benjamini-Hochberg)", value = TRUE),
    hr(),
    section_note("Unit: 51 local municipalities of KwaZulu-Natal (Census 2011).",
                 "Weights are row-standardised. Baseline for the presentation: queen contiguity, α = 0.05, all four predictors.")
  ),

  # ---- Home ----
  nav_panel("Home",
    div(class = "hero",
      h2("Unemployment and Deprivation in Education and Basic Services in KwaZulu-Natal: Associations and Spatial Dependence"),
      div(class = "lbl", "RESEARCH QUESTION"),
      p(em("Is unemployment linked to deprivation in education and basic services across KwaZulu-Natal, and is this relationship spatial?")),
      div(class = "lbl", "AIM"),
      p("To assess the association between unemployment and deprivation in education and basic services across the 51 local municipalities of KwaZulu-Natal, and to test for spatial dependence.", class = "mb-0")
    ),
    layout_column_wrap(width = 1/3, fill = FALSE,
      obj_card("1", "Association", "Examine links between unemployment, schooling and access to water and electricity using correlations and regression.", gold),
      obj_card("2", "Patterns", "Map municipal differences and test spatial clustering using global Moran's I and LISA.", blue),
      obj_card("3", "Spatial effects", "Test whether spatial dependence remains after controlling for deprivation, using OLS, SLM and SEM.", teal)
    ),
    layout_column_wrap(width = 1/4, fill = FALSE,
      value_box("Mean municipal rate (unweighted)", textOutput("vb_mean"), theme = "primary"),
      value_box("Highest", textOutput("vb_max"), textOutput("vb_max_name"), theme = "danger"),
      value_box("Lowest", textOutput("vb_min"), textOutput("vb_min_name"), theme = "success"),
      value_box("Global Moran's I", textOutput("vb_moran"), textOutput("vb_moran_p"), theme = "secondary")
    ),
    card(card_header("Baseline results (queen contiguity, α = 0.05, all four predictors)"),
      tags$ul(
        tags$li("Unemployment is associated with deprivation; the strongest link is piped water (r = −0.67)."),
        tags$li("Unemployment itself shows no significant global clustering (Moran's I = −0.022, p = 0.51), but no schooling and electricity are clustered."),
        tags$li("OLS residuals are spatially dependent (Moran's I = 0.199, p = 0.003); the spatial error model has the lowest AIC (340.4)."),
        tags$li(tags$b("The relationship has a spatial dimension, mainly through spatially correlated unobserved factors.")))),
    layout_columns(col_widths = c(4, 4, 4),
      card(card_header("Unemployment"), p("South Africa's most pressing socio-economic problem. In 2011 the KZN rate ranged from 16% to 57% between municipalities.")),
      card(card_header("Deprivation"), p("Unemployment is one dimension of multidimensional poverty, alongside education and living standards (water, electricity) in Stats SA's SAMPI.")),
      card(card_header("Space"), p("Tobler's first law: nearby places tend to be similar. Classical statistics assumes independence, so we test for spatial dependence."))
    )
  ),

  # ---- Objective 1: Association ----
  nav_panel("1 · Association",
    navset_card_underline(
      nav_panel("Map",
        layout_columns(col_widths = c(8, 4),
          card(full_screen = TRUE, card_header(textOutput("map_title")),
               leafletOutput("map", height = 540),
               card_footer(class = "note", "Hover for the value; click a municipality to highlight it everywhere.",
                           downloadButton("dl_map", "Map (PNG)", class = "btn-sm"))),
          card(full_screen = TRUE, card_header("Ranking of municipalities"),
               plotlyOutput("rank", height = 580))
        )
      ),
      nav_panel("Distribution",
        layout_columns(col_widths = c(6, 6),
          card(card_header("Histogram"), plotlyOutput("hist", height = 320)),
          card(card_header("Box plot"), plotlyOutput("box", height = 320))
        ),
        layout_columns(col_widths = c(6, 6),
          card(card_header("Normal QQ plot"), plotlyOutput("qq", height = 320)),
          card(card_header("Summary statistics (unweighted, across 51 municipalities)"), DTOutput("summ"), uiOutput("skew_text"))
        ),
        card(card_header("Shapiro-Wilk test for normality (descriptive diagnostic)"), uiOutput("shapiro"))
      ),
      nav_panel("Correlation",
        section_note(tags$b("Reading the indicators: "),
                     "higher matric, piped water and electricity percentages mean ", tags$b("less"),
                     " deprivation; a higher no-schooling percentage means ", tags$b("more"), " deprivation."),
        layout_columns(col_widths = c(6, 6),
          card(card_header("Unemployment rate (%)"), leafletOutput("map_u", height = 380),
               card_footer(class = "note", "The right-hand map follows this map's zoom and pan.")),
          card(card_header(selectInput("var2", NULL, choices = vars[-1], selected = "piped_water", width = "100%")),
               leafletOutput("map_v", height = 380))
        ),
        layout_columns(col_widths = c(7, 5),
          card(card_header("Relationship with unemployment"), plotlyOutput("scatter", height = 360), uiOutput("cor_text")),
          card(card_header("Pearson correlation with unemployment (n = 51)"), DTOutput("cor_table"),
               h6("Why correlation is not causation here", class = "obj-title mt-2"),
               tags$ul(
                 tags$li(tags$b("Spatial confounding: "), "shared regional conditions (rurality, distance to towns, former homeland areas) can produce both poor services and high unemployment."),
                 tags$li(tags$b("Non-independence: "), "Pearson p-values assume independent observations. Neighbouring municipalities have similar deprivation, so 51 municipalities carry less information than 51 independent observations, and the p-values are too optimistic."),
                 tags$li(tags$b("Ecological fallacy: "), "an association between municipalities need not hold for individual households.")
               ))
        )
      )
    )
  ),

  # ---- Objective 2: Patterns ----
  nav_panel("2 · Patterns",
    navset_card_underline(
      nav_panel("Spatial weights",
        layout_columns(col_widths = c(7, 5),
          card(card_header("Defining neighbours"),
               tags$ul(
                 tags$li(tags$b("Queen contiguity: "), "neighbours share a border or a corner (main specification)."),
                 tags$li(tags$b("Rook contiguity: "), "neighbours share a border only."),
                 tags$li(tags$b("k-nearest neighbours: "), "the k closest municipality representative points (great-circle distance).")
               ),
               p(tags$b("Row-standardised (style \"W\"): "), "each row of W sums to 1, so the spatial lag Wy is the average of the neighbours' values.")),
          card(card_header("Neighbours per municipality (current weights)"), DTOutput("nb_summary"))
        )
      ),
      nav_panel("Global Moran's I",
        layout_columns(col_widths = c(6, 6),
          card(card_header("Moran scatterplot"), plotlyOutput("moranplot", height = 400),
               card_footer(class = "note", "Slope of the line = Moran's I. Hover over a point to see the municipality.")),
          card(card_header("Hypothesis test (analytical, one-sided)"), uiOutput("moran_rek"))
        ),
        layout_columns(col_widths = c(6, 6),
          card(card_header("Monte Carlo permutation test (999 simulations)"), plotlyOutput("mc_hist", height = 320)),
          card(card_header("Monte Carlo result"), uiOutput("mc_rek"))
        ),
        card(card_header("Global Moran's I for every variable (current weights, one-sided: positive autocorrelation)"), DTOutput("moran_all"),
             section_note("Deprivation indicators can be clustered even when unemployment itself is not."))
      ),
      nav_panel("Local clusters (LISA)",
        layout_columns(col_widths = c(8, 4),
          card(full_screen = TRUE, card_header(textOutput("lisa_title")), leafletOutput("lisa_map", height = 520),
               card_footer(downloadButton("dl_lisa_map", "LISA map (PNG)", class = "btn-sm"))),
          card(card_header("Cluster counts"), DTOutput("lisa_counts"),
               section_note("High-High = hot spot, Low-Low = cold spot, High-Low / Low-High = spatial outliers.",
                            "51 local tests can produce chance findings; Benjamini-Hochberg adjustment (sidebar) controls the false discovery rate. Local clusters are exploratory."))
        ),
        card(card_header("Local Moran's I for every municipality (raw and BH-adjusted p-values)"),
             DTOutput("lisa_table"),
             card_footer(downloadButton("dl_lisa", "LISA results (CSV)", class = "btn-sm")))
      ),
      nav_panel("Hot spots (Getis-Ord Gi*)",
        layout_columns(col_widths = c(8, 4),
          card(full_screen = TRUE, card_header(textOutput("gi_title")), leafletOutput("gi_map", height = 520),
               card_footer(class = "note", "Hover a municipality for its Gi* z-score and p-values.",
                           downloadButton("dl_gi_map", "Gi* map (PNG)", class = "btn-sm"))),
          card(card_header("How to read Gi*"),
               p("Gi* compares the average value in a municipality ", tags$b("and its neighbours"),
                 " with the provincial average. It gives one continuous z-score instead of four quadrant classes:"),
               tags$ul(tags$li(tags$b("Positive z: "), "a concentration of high values (hot spot)."),
                       tags$li(tags$b("Negative z: "), "a concentration of low values (cold spot)."),
                       tags$li(tags$b("Near 0: "), "no local concentration.")),
               p(class = "note", "Confidence classes: |z| ≥ 1.65 (90%), ≥ 1.96 (95%), ≥ 2.58 (99%), two-sided. ",
                 "When Benjamini-Hochberg is ticked in the sidebar, classes use adjusted p-values instead."),
               DTOutput("gi_counts"))
        ),
        card(card_header("Gi* for every municipality"), DTOutput("gi_table"),
             card_footer(downloadButton("dl_gi", "Gi* results (CSV)", class = "btn-sm")))
      )
    )
  ),

  # ---- Objective 3: Spatial effects ----
  nav_panel("3 · Spatial effects",
    uiOutput("scenario_note_reg"),
    card(card_header("Model: unemployment rate explained by"),
         checkboxGroupInput("preds", NULL, choices = predictors, selected = predictors, inline = TRUE)),
    layout_columns(col_widths = c(6, 6),
      card(card_header("OLS coefficients"), DTOutput("ols"), uiOutput("ols_fit")),
      card(card_header("Collinearity (variance inflation factors)"), DTOutput("vif"),
           section_note("A VIF above 5 suggests the coefficient is hard to interpret on its own because that predictor overlaps with the others."))
    ),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Are the OLS residuals spatially autocorrelated?"), uiOutput("resid_rek")),
      card(card_header("OLS residual map"), leafletOutput("resid_map", height = 380))
    ),
    card(card_header("Lagrange Multiplier tests: which spatial model?"), DTOutput("lm_tests"), uiOutput("lm_choice")),
    card(card_header("Model comparison: OLS vs spatial lag (SLM) vs spatial error (SEM)"),
         DTOutput("model_comp"), uiOutput("model_text"),
         card_footer(downloadButton("dl_models", "Model comparison (CSV)", class = "btn-sm"))),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Spatial lag model coefficients"), DTOutput("slm_coef")),
      card(card_header("Spatial error model coefficients"), DTOutput("sem_coef"))
    )
  ),

  # ---- Analysis and conclusion ----
  nav_panel("Analysis & conclusion",
    uiOutput("scenario_note"),
    layout_columns(col_widths = c(4, 4, 4),
      card(card_header("1. Compare variables"),
           tags$ul(tags$li("Paired choropleth maps and scatterplots."),
                   tags$li("Higher unemployment coincides with less water, electricity and schooling."),
                   tags$li("Correlation describes association, not cause."))),
      card(card_header("2. Test clustering"),
           tags$ul(tags$li("Moran's I for each variable; LISA for local patterns."),
                   tags$li("Education and electricity cluster; unemployment shows no significant global clustering."))),
      card(card_header("3. Model spatial effects"),
           tags$ul(tags$li("OLS, then residual Moran's I; compare SLM and SEM with diagnostics and AIC."),
                   tags$li("SEM suggests spatially correlated omitted factors, not proven spillovers.")))
    ),
    div(class = "key", tags$b("Key distinction: "),
        "spatially patterned deprivation ≠ globally clustered unemployment; regression residuals can still be spatially dependent."),
    br(),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Conclusions (current settings)"), uiOutput("conclusions")),
      card(card_header("Limitations"),
           tags$ul(tags$li("Census 2011 is dated; results cannot describe current conditions."),
                   tags$li("Only 51 municipalities reduce statistical power and hide within-municipality inequality (MAUP)."),
                   tags$li("Excluding neighbouring provinces creates boundary effects."),
                   tags$li("Correlated education and service measures complicate coefficient interpretation; associations do not establish causation."),
                   tags$li("Unmeasured factors such as accessibility and local economic structure may influence results.")))
    )
  ),

  # ---- Data ----
  nav_panel("Data",
    card(card_header("Census 2011 indicators by municipality"),
         DTOutput("data_table"),
         card_footer(downloadButton("dl", "Download CSV")))
  ),

  # ---- About ----
  nav_panel("About",
    layout_columns(col_widths = c(7, 5),
      card(card_header("Methods and data"),
        markdown("
**Data.** Census 2011 indicators: Statistics South Africa (2012) *Census 2011 Municipal Fact Sheet*.
Boundaries: Municipal Demarcation Board, 2011 local municipalities (51 in KZN).

**Variables.** Outcome: official unemployment rate (youth unemployment shown on maps).
Deprivation indicators (education and living-standards dimensions of Stats SA's multidimensional poverty index):
% of adults (20+) with matric and with no schooling; % of households with piped water inside the dwelling
and using electricity for lighting. Summary statistics are unweighted means across municipalities, not
population-weighted provincial rates.

**Methods.** Exploratory data analysis and Pearson correlation (objective 1); spatial weights, global Moran's I
(analytical and Monte Carlo, one-sided for positive autocorrelation), LISA and Getis-Ord Gi* hot spots, with optional
Benjamini-Hochberg adjustment (objective 2); OLS with VIFs, Moran's I on the residuals, Lagrange Multiplier tests, and spatial lag /
spatial error models compared with AIC and residual Moran's I (objective 3). Every test is reported in four steps:
hypotheses, statistic and p-value, decision, conclusion.

**Course.** Spatial Statistics (STAT7006A), University of the Witwatersrand

**Group members.**

1. Lintle Lenepa
2. Sheofatso Padima
3. Fanisa Mabunda
4. Congress Mafukele

**Date.** 26 October 2026

**References.** Anselin (1995); Anselin et al. (1996); Benjamini and Hochberg (1995); Bivand and Wong (2018); Getis and Ord (1992); Ord and Getis (1995);
Moran (1950); Statistics South Africa (2012); Tobler (1970).
")),
      card(card_header("Data checks and software"),
           tableOutput("checks"), tableOutput("versions"))
    )
  )
)

# ===================== SERVER =====================
server <- function(input, output, session) {

  nbr   <- reactive(make_nb(input$wtype, input$k))
  lw    <- reactive(nb2listw(nbr(), style = "W", zero.policy = TRUE))
  wlab  <- reactive(if (input$wtype == "knn") paste0("k-nearest neighbours (k = ", input$k, ")") else w_names[[input$wtype]])
  selected_code <- reactive(if (is.null(input$selected_muni)) "" else input$selected_muni)
  for (map_id in c("map", "map_u", "map_v")) {
    local({
      this_map <- map_id
      observeEvent(input[[paste0(this_map, "_shape_click")]], {
        click <- input[[paste0(this_map, "_shape_click")]]
        if (!is.null(click$id)) updateSelectizeInput(session, "selected_muni", selected = as.character(click$id))
      }, ignoreInit = TRUE)
    })
  }
  alpha <- reactive(as.numeric(input$alpha))
  x     <- reactive(muni[[input$var]])
  mt    <- reactive(moran.test(x(), lw(), zero.policy = TRUE))
  is_baseline <- reactive(input$wtype == "queen" && alpha() == 0.05 &&
                            setequal(input$preds, unname(predictors)))

  # ---- Selected municipality: at a glance ----
  output$glance <- renderUI({
    code <- selected_code(); if (!nzchar(code)) return(NULL)
    r <- muni_df[as.character(muni_df$muni_code) == code, ]
    if (nrow(r) == 0) return(NULL)
    tags$table(class = "glance",
      tags$tr(tags$td(tags$b(r$muni_name), colspan = 2)),
      tags$tr(tags$td("Code / district"), tags$td(paste(r$muni_code, "/", r$district))),
      tags$tr(tags$td("Unemployment"), tags$td(paste0(r$unemployment, "%"))),
      tags$tr(tags$td("Youth unemployment"), tags$td(paste0(r$youth_unemp, "%"))),
      tags$tr(tags$td("No schooling"), tags$td(paste0(r$no_school, "%"))),
      tags$tr(tags$td("Matric"), tags$td(paste0(r$matric, "%"))),
      tags$tr(tags$td("Piped water"), tags$td(paste0(r$piped_water, "%"))),
      tags$tr(tags$td("Electricity"), tags$td(paste0(r$electricity, "%"))))
  })

  scenario_html <- reactive({
    if (is_baseline()) return(NULL)
    div(class = "scenario", tags$b("Interactive scenario: "),
        sprintf("results below use %s weights, α = %s and %d predictor(s). The presentation's baseline is queen contiguity, α = 0.05, all four predictors.",
                wlab(), alpha(), length(input$preds)))
  })
  output$scenario_note <- renderUI(scenario_html())
  output$scenario_note_reg <- renderUI(scenario_html())

  # ---- Overview ----
  output$vb_mean <- renderText(sprintf("%.1f%%", mean(x(), na.rm = TRUE)))
  output$vb_max  <- renderText(sprintf("%.1f%%", max(x())))
  output$vb_max_name <- renderText(muni$muni_name[which.max(x())])
  output$vb_min  <- renderText(sprintf("%.1f%%", min(x())))
  output$vb_min_name <- renderText(muni$muni_name[which.min(x())])
  output$vb_moran <- renderText(sprintf("%.3f", mt()$estimate[1]))
  output$vb_moran_p <- renderText({
    p <- mt()$p.value
    paste0("p ", ifelse(p < 0.001, "< 0.001", paste("=", round(p, 3))),
           ifelse(p < alpha(), " · positive autocorrelation", " · no significant positive autocorrelation"))
  })
  output$map_title <- renderText(paste0(lab(input$var), ", Census 2011"))
  output$map <- renderLeaflet(choropleth(input$var, lab(input$var), selected_code()))
  output$dl_map <- downloadHandler(
    filename = function() paste0("map_", input$var, ".png"),
    content = function(file) {
      v <- x(); br <- pretty(v, 5); cl <- cut(v, br, include.lowest = TRUE)
      cols <- hcl.colors(nlevels(cl), "YlOrRd", rev = TRUE)
      png_map(file, cols[as.integer(cl)], paste0(lab(input$var), ", KwaZulu-Natal"),
              levels(cl), cols, "Source: Statistics South Africa, Census 2011; Municipal Demarcation Board 2011 boundaries.")
    })

  output$rank <- renderPlotly({
    d <- muni_df %>% mutate(val = .data[[input$var]]) %>% arrange(val)
    d$highlight <- ifelse(as.character(d$muni_code) == selected_code(), "#142B50", "#E0A526")
    plot_ly(d, x = ~val, y = ~factor(muni_name, levels = muni_name), type = "bar",
            orientation = "h", marker = list(color = d$highlight),
            hovertemplate = "%{y}: %{x}%<extra></extra>") %>%
      layout(xaxis = list(title = paste0(lab(input$var), ", Census 2011")),
             yaxis = list(title = "", tickfont = list(size = 9)), margin = list(l = 10)) %>%
      pconf(paste0("ranking_", input$var))
  })

  # ---- Compare ----
  output$map_u <- renderLeaflet(choropleth("unemployment", "Unemployment", selected_code()))
  output$map_v <- renderLeaflet(choropleth(input$var2, lab(input$var2), selected_code()))
  observeEvent(input$map_u_bounds, {          # one-way sync: right map follows left map
    b <- input$map_u_bounds
    leafletProxy("map_v") %>% fitBounds(b$west, b$south, b$east, b$north)
  }, ignoreInit = TRUE)
  output$scatter <- renderPlotly({
    d <- muni_df %>% mutate(v = .data[[input$var2]])
    fit <- lm(unemployment ~ v, data = d)
    xs <- seq(min(d$v), max(d$v), length.out = 50)
    d$point_color <- ifelse(as.character(d$muni_code) == selected_code(), "#E0A526", "#3B6CF6")
    d$point_size <- ifelse(as.character(d$muni_code) == selected_code(), 16, 9)
    plot_ly(d, x = ~v, y = ~unemployment, type = "scatter", mode = "markers",
            text = ~muni_name, marker = list(size = d$point_size, color = d$point_color, opacity = 0.85),
            hovertemplate = "%{text}<br>x = %{x}%<br>Unemployment = %{y}%<extra></extra>") %>%
      add_lines(x = xs, y = predict(fit, data.frame(v = xs)), inherit = FALSE,
                line = list(color = "#c0392b"), name = "Linear fit", hoverinfo = "skip") %>%
      layout(showlegend = FALSE, xaxis = list(title = paste0(lab(input$var2), ", Census 2011")),
             yaxis = list(title = "Unemployment rate (%)")) %>%
      pconf(paste0("scatter_unemployment_", input$var2))
  })
  output$cor_text <- renderUI({
    ct <- cor.test(muni_df$unemployment, muni_df[[input$var2]])
    r <- unname(ct$estimate)
    strength <- cut(abs(r), c(-Inf, 0.3, 0.6, Inf), labels = c("weak", "moderate", "strong"))
    dir <- ifelse(r > 0, "positive", "negative")
    HTML(sprintf("<h4>r = %.2f</h4><p>A <b>%s %s</b> linear relationship:
                  as %s increases, unemployment tends to %s.</p>
                  <p class='note'>Pearson test: t = %.2f, df = %d, p = %s (assumes independent observations). Correlation is not causation, and it ignores space.</p>",
                 r, strength, dir, tolower(lab(input$var2)),
                 ifelse(r > 0, "increase", "decrease"), ct$statistic, as.integer(ct$parameter), fmt_p(ct$p.value)))
  })

  # ---- Distribution ----
  output$hist <- renderPlotly({
    v <- x(); br <- pretty(v, 12); cl <- cut(v, br, include.lowest = TRUE)
    n <- as.vector(table(cl)); mids <- (head(br, -1) + tail(br, -1)) / 2
    who <- vapply(levels(cl), function(l) paste(sort(muni$muni_name[cl == l]), collapse = ", "), character(1))
    tip <- sprintf("<b>%s–%s%%</b><br>%d municipalities (%.1f%% of 51)<br><i>%s</i>",
                   head(br, -1), tail(br, -1), n, 100 * n / length(v), ifelse(nzchar(who), who, "none"))
    plot_ly(x = mids, y = n, type = "bar", width = diff(br)[1] * 0.95, text = tip, textposition = "none",
            hovertemplate = "%{text}<extra></extra>",
            marker = list(color = "#3B6CF6", line = list(color = "white", width = 1))) %>%
      layout(xaxis = list(title = lab(input$var)), yaxis = list(title = "Number of municipalities"),
             hoverlabel = list(align = "left")) %>% pconf(paste0("histogram_", input$var))
  })
  output$box <- renderPlotly(
    plot_ly(y = x(), type = "box", name = "", boxpoints = "all", jitter = 0.4,
            text = muni$muni_name, marker = list(color = "#3B6CF6"),
            line = list(color = "#3B6CF6"),
            hovertemplate = "<b>%{text}</b><br>%{y}%<extra></extra>") %>%
      layout(yaxis = list(title = lab(input$var)), hovermode = "closest") %>% pconf(paste0("boxplot_", input$var)))
  output$qq <- renderPlotly({
    q <- qqnorm(x(), plot.it = FALSE)
    slope <- diff(quantile(x(), c(0.25, 0.75))) / diff(qnorm(c(0.25, 0.75)))
    int <- quantile(x(), 0.25) - slope * qnorm(0.25)
    plot_ly(x = q$x, y = q$y, type = "scatter", mode = "markers", text = muni$muni_name,
            marker = list(color = "#3B6CF6", size = 8),
            hovertemplate = "<b>%{text}</b><br>Observed: %{y}%<br>Theoretical quantile: %{x:.2f}<extra></extra>") %>%
      add_lines(x = range(q$x), y = int + slope * range(q$x), inherit = FALSE,
                line = list(color = "#c0392b"), hoverinfo = "skip") %>%
      layout(showlegend = FALSE, xaxis = list(title = "Theoretical quantiles"),
             yaxis = list(title = "Sample quantiles")) %>% pconf(paste0("qqplot_", input$var))
  })
  output$summ <- renderDT(dt_small({
    v <- x()
    data.frame(Statistic = c("n", "Mean", "Median", "Standard deviation", "IQR", "Minimum", "Maximum"),
               Value = c(length(v), mean(v), median(v), sd(v), IQR(v), min(v), max(v)))
  }, 2))
  output$skew_text <- renderUI({
    m <- mean(x()); md <- median(x())
    txt <- if (abs(m - md) < 0.5) "approximately symmetric" else if (m > md)
      "positively (right) skewed" else "negatively (left) skewed"
    HTML(sprintf("<p>Mean (%.1f) %s median (%.1f), so the distribution is <b>%s</b>.</p>",
                 m, ifelse(abs(m - md) < 0.5, "≈", ifelse(m > md, "&gt;", "&lt;")), md, txt))
  })
  output$shapiro <- renderUI({
    sw <- shapiro.test(x())
    why <- if (sw$p.value >= alpha()) sprintf(
      "<p class='note'><b>Why we fail to reject:</b> a W this far from 1 would occur about %.0f%% of the time if the data were truly normal, so the departure is too small to count as evidence. Failing to reject does not prove normality; it means we have no evidence against it. The permutation (Monte Carlo) Moran test does not rely on normality.</p>",
      100 * sw$p.value) else ""
    tagList(rek_box(paste0("the ", tolower(lab(input$var)), " is normally distributed."),
            "it is not normally distributed.",
            sprintf("W = %.3f, p = %s.", sw$statistic, fmt_p(sw$p.value)),
            sw$p.value, alpha(),
            "There is sufficient evidence that the variable is not normally distributed.",
            "There is insufficient evidence against normality; this does not establish normality or validate every model assumption."),
            HTML(why))
  })

  # ---- Global Moran's I ----
  output$moranplot <- renderPlotly({
    z <- as.numeric(scale(x())); wz <- lag.listw(lw(), z, zero.policy = TRUE)
    b <- coef(lm(wz ~ z))[2]
    quad <- ifelse(z >= 0 & wz >= 0, "High-High", ifelse(z < 0 & wz < 0, "Low-Low",
            ifelse(z >= 0, "High-Low", "Low-High")))
    plot_ly(x = z, y = wz, type = "scatter", mode = "markers", color = quad,
            colors = lisa_cols[1:4], text = muni$muni_name, marker = list(size = 9),
            hovertemplate = "%{text}<br>z = %{x:.2f}<br>lag = %{y:.2f}<extra></extra>") %>%
      add_lines(x = range(z), y = b * range(z), inherit = FALSE,
                line = list(color = "black"), name = sprintf("slope = I = %.3f", b)) %>%
      layout(xaxis = list(title = paste("Standardised", tolower(lab(input$var))), zeroline = TRUE),
             yaxis = list(title = paste0("Spatially lagged value (", wlab(), ")"), zeroline = TRUE),
             legend = list(orientation = "h", y = -0.2)) %>% pconf(paste0("moran_scatter_", input$var))
  })
  output$moran_rek <- renderUI({
    m <- mt(); est <- m$estimate
    rek_box("there is no spatial autocorrelation in the values.",
            "there is positive spatial autocorrelation (similar values cluster). One-sided test: negative autocorrelation is not tested.",
            sprintf("Moran's I = %.4f, E[I] = −1/(n−1) = %.4f, z = %.3f, p = %s (%s weights).",
                    est[1], est[2], m$statistic, fmt_p(m$p.value), wlab()),
            m$p.value, alpha(),
            paste0("There is sufficient evidence of positive spatial autocorrelation: municipalities with similar ",
                   tolower(lab(input$var)), " are clustered in KwaZulu-Natal."),
            paste0("There is insufficient evidence of positive spatial autocorrelation for ",
                   tolower(lab(input$var)), "; this does not prove spatial randomness."))
  })
  mc <- reactive({ set.seed(2026); moran.mc(x(), lw(), nsim = 999, zero.policy = TRUE) })
  output$mc_hist <- renderPlotly({
    sims <- mc()$res[-length(mc()$res)]; obs <- mc()$statistic
    h <- hist(sims, breaks = 40, plot = FALSE)
    tip <- sprintf("I from %.3f to %.3f<br>%d of 999 permutations (%.1f%%)",
                   head(h$breaks, -1), tail(h$breaks, -1), h$counts, 100 * h$counts / 999)
    ymax <- max(h$counts) * 1.05
    plot_ly(x = h$mids, y = h$counts, type = "bar", width = diff(h$breaks)[1] * 0.95, name = "Simulated I",
            text = tip, textposition = "none", hovertemplate = "%{text}<extra></extra>",
            marker = list(color = "#9fb8c7")) %>%
      add_lines(x = c(obs, obs), y = c(0, ymax), inherit = FALSE, name = "Observed I",
                line = list(color = "#c0392b", width = 3),
                hovertemplate = sprintf("<b>Observed I = %.4f</b><br>pseudo p = %s<extra></extra>", obs, fmt_p(mc()$p.value))) %>%
      layout(xaxis = list(title = "Moran's I under random permutation"),
             yaxis = list(title = "Number of permutations"), legend = list(orientation = "h", y = -0.25)) %>%
      pconf(paste0("moran_permutation_", input$var))
  })

  output$mc_rek <- renderUI({
    m <- mc()
    n_ge <- sum(m$res[-length(m$res)] >= m$statistic)
    rek_box("there is no positive spatial autocorrelation under the permutation null.",
            "there is positive spatial autocorrelation (alternative: greater).",
            sprintf("Observed I = %.4f; %d of 999 permuted values were ≥ the observed I; pseudo p-value = (%d + 1)/(999 + 1) = %s.",
                    m$statistic, n_ge, n_ge, fmt_p(m$p.value)),
            m$p.value, alpha(),
            "The observed clustering is very unlikely under random rearrangement, so there is sufficient evidence of positive spatial autocorrelation.",
            "The observed Moran's I is typical of random rearrangements, so there is insufficient evidence of positive spatial autocorrelation.")
  })

  # ---- LISA ----
  lisa <- reactive({
    li <- localmoran(x(), lw(), zero.policy = TRUE)
    z <- as.numeric(scale(x())); wz <- lag.listw(lw(), z, zero.policy = TRUE)
    p_raw <- li[, 5]
    p <- if (isTRUE(input$lisa_fdr)) p.adjust(p_raw, method = "BH") else p_raw
    quad <- ifelse(z > 0 & wz > 0, "High-High", ifelse(z < 0 & wz < 0, "Low-Low",
            ifelse(z > 0, "High-Low", "Low-High")))          # Moran scatterplot quadrant
    cl <- ifelse(is.na(p) | p >= alpha(), "Not significant", quad)
    data.frame(Municipality = muni$muni_name, District = muni$district,
               Value = x(), Ii = round(li[, 1], 3), Quadrant = quad, `Raw p` = signif(p_raw, 3),
               `BH-adjusted p` = signif(p.adjust(p_raw, method = "BH"), 3),
               Cluster = factor(cl, levels = names(lisa_cols)), check.names = FALSE)
  })
  output$lisa_title <- renderText(paste0("LISA cluster map: ", lab(input$var), " (",
                                         if (isTRUE(input$lisa_fdr)) "BH-adjusted" else "unadjusted",
                                         " p < ", alpha(), "; ", wlab(), ")"))
  output$lisa_map <- renderLeaflet({
    d <- lisa(); pal <- colorFactor(lisa_cols, levels = names(lisa_cols))
    p_show <- if (isTRUE(input$lisa_fdr)) d$`BH-adjusted p` else d$`Raw p`
    base_map() %>%
      addPolygons(fillColor = pal(d$Cluster), fillOpacity = 0.85, color = "white", weight = 0.7,
                  label = lapply(paste0("<b>", d$Municipality, "</b><br>Map class: <b>", d$Cluster, "</b>",
                                        "<br>Quadrant: ", d$Quadrant,
                                        ifelse(d$Cluster == "Not significant",
                                               paste0(" (pattern not significant at p < ", alpha(), ")"), ""),
                                        "<br>Local I = ", d$Ii, "; tested p = ", p_show,
                                        "<br>Raw p = ", d$`Raw p`, "; BH-adjusted p = ", d$`BH-adjusted p`), HTML),
                  popup = popup_html,
                  highlightOptions = highlightOptions(weight = 3, color = "#333", bringToFront = TRUE)) %>%
      addLegend("bottomright", colors = unname(lisa_cols), labels = names(lisa_cols),
                title = paste0("LISA cluster<br>p < ", alpha(), if (isTRUE(input$lisa_fdr)) " (BH)" else " (raw)")) %>%
      add_selection_outline(selected_code())
  })
  output$lisa_counts <- renderDT(dt_small(as.data.frame(table(Cluster = lisa()$Cluster)) %>%
                                      rename(Municipalities = Freq), 2))
  output$lisa_table <- renderDT(
    datatable(lisa(), rownames = FALSE, filter = "top",
              options = list(pageLength = 10, order = list(list(5, "asc")))))
  output$dl_lisa <- downloadHandler(
    filename = function() paste0("lisa_", input$var, "_", input$wtype, ".csv"),
    content = function(file) write.csv(lisa(), file, row.names = FALSE))
  output$dl_lisa_map <- downloadHandler(
    filename = function() paste0("lisa_map_", input$var, ".png"),
    content = function(file) {
      d <- lisa()
      png_map(file, unname(lisa_cols[as.character(d$Cluster)]),
              paste0("LISA clusters: ", lab(input$var)), names(lisa_cols), unname(lisa_cols),
              paste0("Census 2011; ", wlab(), " weights; ",
                     if (isTRUE(input$lisa_fdr)) "BH-adjusted" else "unadjusted", " p < ", alpha(), "."))
    })

  # ---- Getis-Ord Gi* hot spots ----
  gi <- reactive({
    lw_star <- nb2listw(include.self(nbr()), style = "W", zero.policy = TRUE)   # Gi* includes the municipality itself
    z <- as.numeric(localG(x(), lw_star, zero.policy = TRUE))
    p_raw <- 2 * pnorm(-abs(z))
    p_bh <- p.adjust(p_raw, method = "BH")
    p <- if (isTRUE(input$lisa_fdr)) p_bh else p_raw
    lvl <- ifelse(p < 0.01, "99%", ifelse(p < 0.05, "95%", ifelse(p < 0.10, "90%", NA)))
    cls <- ifelse(is.na(lvl), "Not significant", paste0(ifelse(z > 0, "Hot spot (", "Cold spot ("), lvl, ")"))
    data.frame(Code = muni$muni_code, Municipality = muni$muni_name, District = muni$district,
               Value = x(), `Gi* z` = round(z, 3), `Raw p` = signif(p_raw, 3), `BH-adjusted p` = signif(p_bh, 3),
               Class = factor(cls, levels = names(gi_cols)), check.names = FALSE)
  })
  output$gi_title <- renderText(paste0("Getis-Ord Gi* hot spots: ", lab(input$var), " (",
                                       if (isTRUE(input$lisa_fdr)) "BH-adjusted" else "unadjusted",
                                       " p; ", wlab(), ")"))
  output$gi_map <- renderLeaflet({
    d <- gi(); pal <- colorFactor(gi_cols, levels = names(gi_cols))
    base_map() %>%
      addPolygons(fillColor = pal(d$Class), fillOpacity = 0.85, color = "#9aa5b1", weight = 0.7,
                  label = lapply(paste0("<b>", d$Municipality, "</b><br>Class: <b>", d$Class, "</b>",
                                        "<br>Gi* z = ", d$`Gi* z`,
                                        ifelse(d$`Gi* z` > 0, " (high values around it)", " (low values around it)"),
                                        "<br>Raw p = ", d$`Raw p`, "; BH-adjusted p = ", d$`BH-adjusted p`,
                                        "<br>", lab(input$var), ": ", d$Value, "%"), HTML),
                  popup = popup_html,
                  highlightOptions = highlightOptions(weight = 3, color = "#333", bringToFront = TRUE)) %>%
      addLegend("bottomright", colors = unname(gi_cols), labels = names(gi_cols),
                title = paste0("Gi* class", if (isTRUE(input$lisa_fdr)) " (BH)" else " (raw p)")) %>%
      add_selection_outline(selected_code())
  })
  output$gi_counts <- renderDT(dt_small(as.data.frame(table(Class = gi()$Class)) %>% rename(Municipalities = Freq), 0))
  output$gi_table <- renderDT(datatable(gi(), rownames = FALSE, filter = "top",
                                        options = list(pageLength = 10, order = list(list(4, "desc")))))
  output$dl_gi <- downloadHandler(
    filename = function() paste0("gi_star_", input$var, "_", input$wtype, ".csv"),
    content = function(file) write.csv(gi(), file, row.names = FALSE))
  output$dl_gi_map <- downloadHandler(
    filename = function() paste0("gi_star_map_", input$var, ".png"),
    content = function(file) {
      d <- gi()
      png_map(file, unname(gi_cols[as.character(d$Class)]),
              paste0("Getis-Ord Gi* hot spots: ", lab(input$var)), names(gi_cols), unname(gi_cols),
              paste0("Census 2011; ", wlab(), " weights (self included); ",
                     if (isTRUE(input$lisa_fdr)) "BH-adjusted" else "unadjusted", " two-sided p."))
    })

  # ---- Regression ----
  form <- reactive({
    req(length(input$preds) > 0)
    reformulate(input$preds, response = "unemployment")
  })
  ols <- reactive(lm(form(), data = muni_df))
  output$ols <- renderDT(dt_small({
    s <- summary(ols())$coefficients
    data.frame(Term = rownames(s), Estimate = s[, 1], `Std. error` = s[, 2],
               t = s[, 3], `p-value` = fmt_p(s[, 4]), check.names = FALSE)
  }, 3))
  output$vif <- renderDT(dt_small(vif_table(input$preds, muni_df), 2))
  output$ols_fit <- renderUI({
    s <- summary(ols())
    HTML(sprintf("<p>R² = %.3f, adjusted R² = %.3f, AIC = %.1f. Coefficient signs describe conditional associations, not causal effects.
                  Check collinearity (VIF) and model assumptions before interpreting individual coefficients.</p>",
                 s$r.squared, s$adj.r.squared, AIC(ols())))
  })
  resid_mt <- reactive(lm.morantest(ols(), lw(), zero.policy = TRUE))
  output$resid_rek <- renderUI({
    m <- resid_mt()
    rek_box("the OLS residuals are not spatially autocorrelated.",
            "the OLS residuals are positively spatially autocorrelated.",
            sprintf("Moran's I (residuals) = %.4f, z = %.3f, p = %s (%s weights).",
                    m$estimate[1], m$statistic, fmt_p(m$p.value), wlab()),
            m$p.value, alpha(),
            "The residuals are spatially clustered, so the OLS assumption of independent errors is violated. A spatial regression model is needed.",
            "There is insufficient evidence of positive residual spatial autocorrelation; check other diagnostics before accepting OLS.")
  })
  output$resid_map <- renderLeaflet({
    r <- residuals(ols())
    pal <- colorBin("RdBu", r, bins = 5, reverse = TRUE)
    base_map() %>%
      addPolygons(fillColor = pal(r), fillOpacity = 0.85, color = "white", weight = 0.7,
                  label = paste0(muni$muni_name, ": residual ", round(r, 1))) %>%
      addLegend("bottomright", pal = pal, values = r, title = "OLS residual<br>(percentage points)") %>%
      add_selection_outline(selected_code())
  })
  lmt <- reactive({
    if (exists("lm.RStests", envir = asNamespace("spdep"))) {
      spdep::lm.RStests(ols(), lw(), test = "all", zero.policy = TRUE)
    } else {
      spdep::lm.LMtests(ols(), lw(), test = "all", zero.policy = TRUE)
    }
  })
  output$lm_tests <- renderDT(dt_small({
    t <- lmt()
    data.frame(Test = names(t),
               Statistic = sapply(t, function(e) unname(e$statistic)),
               `p-value` = sapply(t, function(e) fmt_p(e$p.value)), check.names = FALSE)
  }, 3))
  output$lm_choice <- renderUI({
    t <- lmt(); p <- sapply(t, function(e) e$p.value); a <- alpha()
    pick <- function(nm) {             # test names differ between spdep versions
      z <- p[names(p) %in% nm]
      if (length(z) == 0L) return(NA_real_)
      unname(z[1])
    }
    err  <- pick(c("LMerr", "RSerr"));     lag  <- pick(c("LMlag", "RSlag"))
    rerr <- pick(c("RLMerr", "adjRSerr")); rlag <- pick(c("RLMlag", "adjRSlag"))
    msg <- if (anyNA(c(err, lag, rerr, rlag))) {
      "LM test names vary across spdep versions; inspect the table and compare model diagnostics and AIC."
    } else if (err >= a && lag >= a) {
      "Neither standard LM test is significant at the selected α. This does not prove that OLS is adequate."
    } else if (err < a && lag >= a) {
      "Only the standard error test is significant, so the SEM is the candidate model. Confirm with AIC and residual diagnostics below."
    } else if (lag < a && err >= a) {
      "Only the standard lag test is significant, so the SLM is the candidate model. Confirm with AIC and residual diagnostics below."
    } else if (rerr < a && rlag < a) {
      "Both standard and both robust tests are significant: the evidence is ambiguous. Compare SEM and SLM using AIC, residual diagnostics and theory."
    } else if (rerr < a) {
      "Only robust error is significant: SEM is supported, subject to fit and residual diagnostics."
    } else if (rlag < a) {
      "Only robust lag is significant: SLM is supported, subject to fit and residual diagnostics."
    } else {
      "Neither robust test is significant. Compare candidate models using theory and diagnostics."
    }
    note <- if (!anyNA(c(rerr, rlag)) && rerr < a && rlag < a)
      sprintf("<p class='note'>Robust LM-error p = %s and robust LM-lag p = %s are both reported; the choice is not made on the smaller p-value alone.</p>",
              fmt_p(rerr), fmt_p(rlag)) else ""
    HTML(paste0("<p><b>Decision rule (Anselin et al., 1996):</b> start with the standard LM-error and LM-lag tests; use the robust versions only when both are significant.</p><p>", msg, "</p>", note))
  })
  models <- reactive({
    d <- muni_df
    list(OLS = ols(),
         SLM = lagsarlm(form(), data = d, listw = lw(), zero.policy = TRUE),
         SEM = errorsarlm(form(), data = d, listw = lw(), zero.policy = TRUE))
  })
  model_table <- reactive({
    m <- models()
    rmi <- lapply(m, function(z) moran.test(residuals(z), lw(), zero.policy = TRUE))
    data.frame(Model = c("OLS", "Spatial lag (SLM)", "Spatial error (SEM)"),
               `Spatial parameter` = c("-", sprintf("ρ = %.3f (LR p %s)", m$SLM$rho,
                                                    fmt_p(summary(m$SLM)$LR1$p.value)),
                                       sprintf("λ = %.3f (LR p %s)", m$SEM$lambda,
                                               fmt_p(summary(m$SEM)$LR1$p.value))),
               `Log-likelihood` = round(sapply(m, function(z) as.numeric(logLik(z))), 2),
               AIC = round(sapply(m, AIC), 2),
               `Residual Moran's I` = sapply(rmi, function(r) sprintf("%.3f (p %s)", r$estimate[1], fmt_p(r$p.value))),
               check.names = FALSE, row.names = NULL)
  })
  output$model_comp <- renderDT(dt_small(model_table(), 2))
  output$dl_models <- downloadHandler(
    filename = function() paste0("model_comparison_", input$wtype, ".csv"),
    content = function(file) write.csv(model_table(), file, row.names = FALSE))
  coef_tab <- function(mod) {
    cf <- summary(mod)$Coef
    data.frame(Term = rownames(cf), Estimate = cf[, 1], `Std. error` = cf[, 2], z = cf[, 3],
               `p-value` = fmt_p(cf[, 4]), check.names = FALSE)
  }
  output$slm_coef <- renderDT(dt_small(coef_tab(models()$SLM), 3))
  output$sem_coef <- renderDT(dt_small(coef_tab(models()$SEM), 3))
  output$model_text <- renderUI({
    a <- sapply(models(), AIC); best <- names(which.min(a))
    full <- c(OLS = "OLS", SLM = "spatial lag model", SEM = "spatial error model")[best]
    HTML(sprintf("<p>The model with the lowest AIC is the <b>%s</b> (AIC = %.1f): it has the best relative fit among the three candidates, not proof of causation.
                  The residual Moran's I column checks whether each model has removed the spatial dependence; AIC alone does not show this.
                  A spatial-lag parameter reflects dependence between neighbouring outcomes but does not by itself establish causal spillovers.
                  A spatial-error parameter indicates spatially correlated unobserved influences.</p>",
                 full, min(a)))
  })

  # ---- Correlation table, weights summary, Moran's I for all, conclusions ----
  output$cor_table <- renderDT(dt_small({
    do.call(rbind, lapply(predictors, function(v) {
      ct <- cor.test(muni_df$unemployment, muni_df[[v]]); r <- unname(ct$estimate)
      data.frame(Variable = lab(v), r = sprintf("%+.2f", r), `p-value` = fmt_p(ct$p.value),
                 Strength = paste(as.character(cut(abs(r), c(-Inf, 0.3, 0.6, Inf), labels = c("weak", "moderate", "strong"))),
                                  ifelse(r > 0, "positive", "negative")), check.names = FALSE)
    }))
  }, 2))
  output$nb_summary <- renderDT(dt_small({
    k <- spdep::card(lw()$neighbours)
    data.frame(Statistic = c("Weights", "Mean", "Minimum", "Maximum", "Municipalities with no neighbours"),
               Value = c(wlab(), round(mean(k), 2), min(k), max(k), sum(k == 0)))
  }, 2))
  output$moran_all <- renderDT(dt_small({
    do.call(rbind, lapply(vars, function(v) {
      m <- moran.test(muni_df[[v]], lw(), zero.policy = TRUE)
      data.frame(Variable = lab(v), `Moran's I` = sprintf("%.3f", m$estimate[1]), `p-value` = fmt_p(m$p.value),
                 `Positive autocorrelation?` = ifelse(m$p.value < alpha(), "Yes", "No"), check.names = FALSE)
    }))
  }, 2))
  output$conclusions <- renderUI({
    r_w <- cor(muni_df$unemployment, muni_df$piped_water)
    mu <- moran.test(muni_df$unemployment, lw(), zero.policy = TRUE)
    rm <- resid_mt(); a <- sapply(models(), AIC); best <- names(which.min(a))
    best_full <- c(OLS = "OLS", SLM = "spatial lag model", SEM = "spatial error model")[best]
    HTML(paste0("<ul>",
      sprintf("<li><b>Objective 1:</b> unemployment is associated with deprivation; the strongest link is piped water (r = %.2f).</li>", r_w),
      sprintf("<li><b>Objective 2:</b> %s of positive global clustering of unemployment (I = %.3f, p = %s).</li>",
              ifelse(mu$p.value < alpha(), "there is evidence", "there is insufficient evidence"), mu$estimate[1], fmt_p(mu$p.value)),
      sprintf("<li><b>Objective 3:</b> OLS residuals %s spatially dependent (I = %.3f, p = %s); the best model by AIC is the %s (AIC = %.1f).</li>",
              ifelse(rm$p.value < alpha(), "are", "are not detectably"), rm$estimate[1], fmt_p(rm$p.value), best_full, min(a)),
      "<li><b>Interpretation:</b> the evidence indicates spatially structured unexplained variation; it does not by itself prove that deprivation causes unemployment.</li></ul>",
      sprintf("<p class='note'>Current settings: %s weights, α = %s, %d predictor(s).</p>", wlab(), alpha(), length(input$preds))))
  })

  # ---- About: checks and versions ----
  output$checks <- renderTable(data_checks)
  output$versions <- renderTable({
    pk <- c("shiny", "sf", "spdep", "spatialreg", "leaflet", "plotly", "DT", "bslib")
    data.frame(Software = c("R", pk),
               Version = c(paste(R.version$major, R.version$minor, sep = "."),
                           sapply(pk, function(p) as.character(utils::packageVersion(p)))))
  })

  # ---- Data ----
  tidy <- muni_df %>% select(Code = muni_code, Municipality = muni_name, District = district,
                             all_of(setNames(vars, names(vars))))
  output$data_table <- renderDT(datatable(tidy, rownames = FALSE, filter = "top",
                                          options = list(pageLength = 15, scrollX = TRUE)))
  output$dl <- downloadHandler("kzn_census2011_unemployment.csv",
                               function(file) write.csv(tidy, file, row.names = FALSE))
}

shinyApp(ui, server)
