#' MMF 2026 Transmission Network Visualisation
#' MMF Exercise
#' MMED 2026, ICI3D Programme
#'
#' This script reads the MMF outbreak data from Google Sheets, cleans it up, and draws an interactive 
#' transmission network using visNetwork.
#'
#' By the end of this script you will have an interactive network where:
#'   Each node is a person in the outbreak
#'   Each arrow shows who infected whom
#'   Node colour tells you whether a person was an index case, symptomatic, or asymptomatic
#'   Node size reflects how many people that person went on to infect
#'
#'

# Loading packages ----------------------------------------------------------------------------

#' We use pacman to load all packages at once. If a package is not installed, pacman installs it 
#' automatically before loading. Each package has a specific job:
#'   dplyr and tidyr: all data manipulation (filtering, selecting, joining)
#'   ggplot2: plotting (not used for the network itself but useful for other plots)
#'   googlesheets4: reading data directly from Google Sheets
#'   lubridate: handling date and time columns cleanly
#'   janitor: cleaning messy column names (removes spaces, lowercases everything)
#'   visNetwork: drawing the interactive transmission network

pacman::p_load(
  dplyr,
  tidyr,
  ggplot2,
  googlesheets4,
  lubridate,
  janitor,
  visNetwork
)

# Importing the data --------------------------------------------------------------------------

#' We authenticate with Google so we can read from Google Sheets. The first time you run this it will
#' open a browser window asking you to log in with your Google account. After that it saves a token 
#' locally so you do not have to log in again every time.
gs4_auth()

#' We read the main exposure data from the "Clean" sheet. This sheet has one row per exposure event, 
#' telling us who infected whom, when the exposure happened, and whether the infected person had 
#' symptoms.
mmf_original <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1ZbxGAdlXW6fjIiKlkj4FiEqRHETT77nzgaTrsOz5J7w",
  sheet = "Clean"
)

#' We read the symptom visit data from the "DoctorVisits" sheet and skip the first row because it is
#' a merged header row that R cannot parse properly. The real column names are on the second row.
symp_original <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1ZbxGAdlXW6fjIiKlkj4FiEqRHETT77nzgaTrsOz5J7w",
  sheet = "DoctorVisits",
  skip = 1
)
#' Alternatively, you can download the data directly and load it from your directory. 
#' If it is already downloaded and in your directory, change the above links to your
#' data directory and load data as any other data. 

# Cleaning the exposure data ------------------------------------------------------------------

#' We clean and simplify the main exposure table. clean_names() converts all column names to lowercase
#' with underscores, so "Person Infected" becomes "person_infected" and so on, then paste date and 
#' time together into one exposure timestamp column because visNetwork will show this as the edge label.
#' We then select and rename only the columns we actually need.
exp <- mmf_original |>
  clean_names() |>
  mutate(expos.time = paste(date, time_sheet)) |>
  select(
    from = infected_by,
    to = person_infected,
    expos.time,
    infection_number,
    symp = symtpomatic
  )

# Building the edges table --------------------------------------------------------------------

#' In network language an "edge" is a connection between two nodes. Here each edge is one transmission 
#' event: person A infected person B. Index cases have no known source so their "from" value is NA.
#' We remove those rows because an edge needs both a from and a to. The arrows column tells visNetwork
#' to draw a directed arrow pointing from the infector to the infected person.
#' The title column is what appears when you hover over an edge.
edges <- exp |>
  filter(!is.na(from)) |>
  mutate(arrows = "to") |>
  rename(title = expos.time)

# Building the nodes table --------------------------------------------------------------------

#' In network language a "node" is a person. We collect every unique person who appears anywhere in
#' the network, either as someone who was infected (the "to" column) or as someone who did the infecting 
#' (the "from" column, excluding NAs).
all_names <- unique(c(exp$to, exp$from[!is.na(exp$from)]))

nodes <- tibble(name = all_names)

#' We then work out who the index cases are. An index case is someone who appears in the "to" column 
#' but whose corresponding "from" value is NA, meaning we do not know who infected them.
index_cases <- exp |>
  filter(is.na(from)) |>
  pull(to)

#' We pull the symptom status for each person from the exposure table. Each person appears as "to" at 
#' least once, and their symp value tells us whether they were symptomatic (1) or asymptomatic (0).
#' We use distinct() to keep only one row per person in case they appear multiple times (for example 
#' if they were exposed more than once).
symp_lookup <- exp |>
  select(name = to, symp) |>
  distinct()

#' We join the symptom information onto the nodes table. left_join keeps everyone in nodes even if 
#' they have no symp record.
nodes <- nodes |>
  left_join(symp_lookup, by = "name")

#' We assign each node a display group based on their status. The group column is what visNetwork uses 
#' to colour nodes. We check index case status first because an index case who was also
#' symptomatic should show up as an index case, not as symptomatic.
nodes <- nodes |>
  mutate(
    group = case_when(
      name %in% index_cases & symp == 1  ~ "Index case (symptomatic)",
      name %in% index_cases & symp == 0  ~ "Index case (asymptomatic)",
      name %in% index_cases ~ "Index case (asymptomatic)",
      symp == 1 ~ "Symptomatic",
      symp == 0 ~ "Asymptomatic",
      is.na(symp) ~ "Asymptomatic"
    )
  )

#' visNetwork requires a column called "id" to uniquely identify each node and a column called "label" 
#' for the text displayed on the node itself. We use the person's name for both.
nodes <- nodes |>
  mutate(
    id = name,
    label = name
  )

#' We calculate the out-degree for each person, which is simply how many people they went on to infect.
#' We count how many times each name appears in the "from" column of the edges table. People who did
#' not infect anyone will be missing from this count, so we fill their N with 0 after joining. We add 
#' 1 to N so that even people who infected nobody still have a visible node size (value = 1). 
#' People who infected many others will have larger nodes, which makes superspreaders immediately visible.
out_degree <- edges |>
  group_by(from) |>
  summarise(n = n()) |>
  rename(name = from)

nodes <- nodes |>
  left_join(out_degree, by = "name") |>
  mutate(
    n = replace_na(n, 0),
    value = n + 1
  )

# Drawing the network -------------------------------------------------------------------------

#' We now pass both tables to visNetwork to draw the interactive network. The nodes table defines who
#' is in the network and how they look. The edges table defines who infected whom. width and height set 
#' the display size of the network in the viewer.
visNetwork(
  nodes,
  edges,
  width = "100%",
  height = "700px",
  main = list(text = "MMF 2026 Transmission Network", style = "font-family:Georgia; font-size:18px; font-weight:bold")
) |>
  #' We style the edges (the arrows between nodes).
  #' color sets the default arrow colour and the colour when selected.
  #' smooth makes the arrows curve slightly so overlapping arrows are
  #' easier to distinguish from one another.
  #' font controls the size and position of the edge label text.
  visEdges(
    arrows = "to",
    color = list(color = "#aaaaaa", highlight = "#e63946"),
    smooth = list(type = "curvedCW", roundness = 0.2),
    font = list(size = 10, align = "middle")
  ) |>
  #' We define the colour for each group we assigned in the nodes table.
  #' background is the fill colour of the node circle.
  #' border is the outline colour.
  #' Red for index cases so they stand out immediately.
  #' Blue for symptomatic cases, orange for asymptomatic, grey for unknown.
  visGroups(
    groupname = "Index case (symptomatic)",
    color = list(background = "#e63946", border = "#9b1a25"),
    shape = "dot"
  ) |>
  visGroups(
    groupname = "Index case (asymptomatic)",
    color = list(background = "#7D3C98", border = "#4a235a"),
    shape = "dot"
  ) |>
  visGroups(
    groupname = "Symptomatic",
    color = list(background = "#2196F3", border = "#0d47a1"),
    shape = "dot"
  ) |>
  visGroups(
    groupname = "Asymptomatic",
    color = list(background = "#FF9800", border = "#e65100"),
    shape = "dot"
  ) |>
  visGroups(
    groupname = "Unknown",
    color = list(background = "#cccccc", border = "#888888"),
    shape = "dot"
  ) |>
  
  #' We add a legend on the right side so the colour coding is clear
  #' to anyone looking at the network for the first time.
  visLegend(position = "right", main = "Node type") |> 
  
  #' We add interactive options.
  #' highlightNearest highlights a node and all its direct and second
  #' degree contacts when you click on it, which is very useful for
  #' tracing transmission chains.
  #' nodesIdSelection adds a dropdown menu so you can search for and
  #' jump to a specific person by name.
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 2, hover = TRUE),
    nodesIdSelection = TRUE
  ) |> 
  #' We fix the random seed so the layout looks the same every time
  #' the script is run. Without this the nodes shuffle into a different
  #' arrangement on each run, which makes it harder to compare across runs.
  visLayout(randomSeed = 152)

#' To save this interative plot, click on Export in the Viewer panel and then click as
#' 'Save as Web Page'.

