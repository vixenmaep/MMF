library(dplyr)
library(readr)
library(ggplot2)
# Make csv files for time windows and hourly
# Visualise Empire v Aims day v night stuff
raw_data <- read_csv("MMF-Final+Locations+Doctors.csv")

# For hours
cleaned_hours <- raw_data %>%
  filter(!is.na(Time) & !is.na(Date)) %>%
  mutate(
    Date = trimws(as.character(Date)),
    Time = trimws(as.character(Time))
  ) %>%
  mutate(Date = ifelse(Date == "18/06", "18/06/2026", Date)) %>%
  mutate(DateTime_Text = paste(Date, Time)) %>%
  mutate(DateTime = as.POSIXct(DateTime_Text, format = "%d/%m/%Y %H:%M")) %>%
  filter(!is.na(DateTime)) %>%
  mutate(Hour_Bin = cut(DateTime, breaks = "hour")) %>%
  mutate(Hour_Bin = as.POSIXct(Hour_Bin))

start_time <- min(cleaned_hours$Hour_Bin, na.rm = TRUE)
end_time   <- max(cleaned_hours$Hour_Bin, na.rm = TRUE)
master_timeline <- data.frame(
  Hour_Bin = seq(from = start_time, to = end_time, by = "hour")
)

hourly_counts <- cleaned_hours %>%
  group_by(Hour_Bin) %>%
  summarise(Count = n(), .groups = 'drop')

full_timeline_data <- master_timeline %>%
  left_join(hourly_counts, by = "Hour_Bin") %>%
  mutate(Count = ifelse(is.na(Count), 0, Count)) %>%
  arrange(Hour_Bin)

ggplot(full_timeline_data, aes(x = Hour_Bin, y = Count)) +
  geom_line(color = "red", size = 1) + 
  geom_area(fill = "tomato", alpha = 0.15) +
  geom_point(data = filter(full_timeline_data, Count > 0), 
             color = "#c0392b", size = 1.5) +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b %d, %Y") +
  theme_minimal(base_size = 13) +
  labs(
    title = "New Infections Per Hour Across All Days",
    subtitle = "Chronological timeline of hourly incident spikes (June 15 - June 19, 2026)",
    x = "Timeline",
    y = "Infection Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b=8)),
    plot.subtitle = element_text(color = "gray40", margin = margin(b=12)),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 15), 
    axis.text = element_text(size = 13)
  )
ggsave("Name of Picture.png", plot label, width = 8, height = 6, dpi = 300)


# For Intervals < 12:30 < 18:30 < 00:00
cleaned_intervals <- raw_data %>%
  filter(!is.na(Time) & !is.na(Date)) %>%
  mutate(
    Date = trimws(as.character(Date)),
    Time = trimws(as.character(Time))
  ) %>%
  mutate(Date = ifelse(Date == "18/06", "18/06/2026", Date)) %>%
  mutate(DateTime_Text = paste(Date, Time)) %>%
  mutate(DateTime = as.POSIXct(DateTime_Text, format = "%d/%m/%Y %H:%M")) %>%
  filter(!is.na(DateTime)) %>%
  
  mutate(
    Hour = as.numeric(format(DateTime, "%H")),
    Minute = as.numeric(format(DateTime, "%M")),
    Time_Minutes = Hour * 60 + Minute
  ) %>%
  mutate(
    Interval_Label = case_when(
      Time_Minutes <= (12 * 60 + 30) ~ "Morning",
      Time_Minutes <= (18 * 60 + 30) ~ "Afternoon",
      TRUE                           ~ "Night"
    )
  ) %>%
  mutate(Date_Label = format(DateTime, "%Y-%m-%d"))

unique_days <- unique(cleaned_intervals$Date_Label)
interval_order <- c("Morning", "Afternoon", "Night")

master_grid <- expand.grid(
  Date_Label = unique_days,
  Interval_Label = interval_order,
  stringsAsFactors = FALSE
) %>%
  arrange(Date_Label, match(Interval_Label, interval_order)) %>%
  mutate(Timeline_Step = paste(Date_Label, Interval_Label))

interval_counts <- cleaned_intervals %>%
  group_by(Date_Label, Interval_Label) %>%
  summarise(Count = n(), .groups = 'drop')

full_timeline_data <- master_grid %>%
  left_join(interval_counts, by = c("Date_Label", "Interval_Label")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count)) %>%
  mutate(Timeline_Step = factor(Timeline_Step, levels = unique(Timeline_Step)))

print(full_timeline_data)

plot <- ggplot(full_timeline_data, aes(x = Timeline_Step, y = Count, group = 1)) +
  # Add a background area shading to visualize the wave shifts
  geom_area(fill = "red", alpha = 0.15) +
  # Continuous trend line across shifts
  geom_line(color = "tomato", size = 1.1) +
  # Highlight the interval endpoints
  geom_point(color = "#c0392b", size = 2.5) +
  
  theme_minimal(base_size = 12) +
  labs(
    title = "Infection Counts Split by Time Windows",
    subtitle = "Window tracking up to 12:30, till 18:30, and till Midnight",
    x = "Daily Time Windows",
    y = "Infection Incident Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b=6)),
    plot.subtitle = element_text(color = "gray40", margin = margin(b=12)),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 15), 
    axis.text = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

ggsave("Data - Time Windows.png", plot, width = 8, height = 6, dpi = 300)

# 1. SET YOUR TARGET X-AXIS POSITIONS HERE
# Replace 3 and 6 with the actual index positions of June 17 Night and June 18 Night
pos_17_night <- 9 
pos_18_night <- 12

max_count <- max(full_timeline_data$Count, na.rm = TRUE)

plot <- ggplot(full_timeline_data, aes(x = Timeline_Step, y = Count, group = 1)) +
  # Add a background area shading to visualize the wave shifts
  geom_area(fill = "red", alpha = 0.15) +
  # Continuous trend line across shifts
  geom_line(color = "tomato", linewidth = 1.1) +
  # Highlight the interval endpoints
  geom_point(color = "#c0392b", size = 2.5) +
  
  # --- Straight Vertical Arrow for Night of 17th June ---
  annotate("segment", 
           x = pos_17_night, y = max_count * 0.85,     # Start directly above the node
           xend = pos_17_night, yend = max_count * 0.70, # Point straight down
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "black") +
  
  # --- Straight Vertical Arrow for Night of 18th June ---
  annotate("segment", 
           x = pos_18_night, y = max_count * 0.85,     # Start directly above the node
           xend = pos_18_night, yend = max_count * 0.70, # Point straight down
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "black") +
  
  # --- Theme & Formatting ---
  theme_minimal(base_size = 12) +
  labs(
    title = "Infection Counts Split by Time Windows",
    subtitle = "Window tracking up to 12:30, till 18:30, and till Midnight",
    x = "Daily Time Windows",
    y = "Infection Incident Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b=6)),
    plot.subtitle = element_text(color = "gray40", margin = margin(b=12)),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 15), 
    axis.text = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

ggsave("Data - Time Windows.png", plot, width = 8, height = 6, dpi = 300)

