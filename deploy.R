# deploy.R - publish the dashboard to shinyapps.io (click "Source")
# Updates https://lintlelenepa.shinyapps.io/kzn-unemployment/
# Run "Run App" on app.R first to check it works on your laptop.
rsconnect::deployApp(appFiles = c("app.R", "data/municipalities.rds"),
                     appName = "kzn-unemployment",
                     forceUpdate = TRUE)
