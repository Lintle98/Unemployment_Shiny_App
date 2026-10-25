# 02_push_to_github.R - put the dashboard on GitHub (run step by step)
# Select ONE step's lines and press Ctrl+Enter. Do not Source the whole file.

# STEP 1 - packages
install.packages(c("usethis", "gitcreds"))

# STEP 2 - tell Git who you are (once per laptop)
usethis::use_git_config(user.name = "Lintle Lenepa",
                        user.email = "lintlelenepa@gmail.com")

# STEP 3 - make this folder a Git repository and make the first commit
#          Answer the prompts: choose the option that means "Yes" to commit,
#          then "Yes" to restart RStudio.
usethis::use_git()

# ---- After RStudio restarts, reopen this file and continue ----

# STEP 4 - create a GitHub token (opens github.com in your browser)
#          Sign in yourself, keep the default boxes ticked, click
#          "Generate token" and COPY the token.
usethis::create_github_token()

# STEP 5 - store the token: paste it when R asks (you type this, not Claude)
gitcreds::gitcreds_set()

# STEP 6 - create the repository on GitHub and push everything
usethis::use_github(private = FALSE)   # use TRUE for a private repo

# Later changes: use the Git tab (top-right) -> tick files -> Commit -> Push
