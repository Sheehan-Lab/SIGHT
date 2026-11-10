# ==============================================================================
# Deploy SIGHT Shiny App to shinyapps.io (https://jkonarkowski.shinyapps.io/sightapp/)
# ==============================================================================
library(rsconnect)

# If your shinyapps.io account token is not yet registered on this machine,
# uncomment and run the line below with your token & secret from:
# https://www.shinyapps.io/admin/#/tokens
#
# rsconnect::setAccountInfo(name = "jkonarkowski", token = "<YOUR_TOKEN>", secret = "<YOUR_SECRET>")

if (nrow(rsconnect::accounts()) == 0) {
  stop("No shinyapps.io account registered yet. Run rsconnect::setAccountInfo(name='jkonarkowski', token='...', secret='...') first.")
}

rsconnect::deployApp(
  appDir = ".",
  appName = "sightapp",
  account = "jkonarkowski",
  forceUpdate = TRUE
)
