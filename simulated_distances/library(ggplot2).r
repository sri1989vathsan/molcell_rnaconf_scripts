library(ggplot2)
library(dplyr)
library(ggplot2)

# Sample data
set.seed(123)
df <- data.frame(
    x = rnorm(100),
    y = rnorm(100),
    group = rep(c("A", "B", "C"), length.out = 100)
)

# Calculate some summary statistics for positioning the text
position_data <- df %>%
    group_by(group) %>%
    summarise(
        median_x = median(x),
        median_y = median(y),
        min_x = min(x),
        max_y = max(y),
        .groups = 'drop'
    ) %>%
    # Create a label for each group
    mutate(label = paste0("Group: ", group, "\nMedian X: ", round(median_x, 2), "\nMedian Y: ", round(median_y, 2)))


# Custom positioning for each group
position_data <- position_data %>%
    mutate(
        text_x = case_when(
            group == "A" ~ min_x - 0.5,  # Position text left of min x for group A
            group == "B" ~ median_x,      # Centered on median x for group B
            group == "C" ~ median_x + 0.5 # Right of median x for group C
        ),
        text_y = case_when(
            group == "A" ~ max_y + 0.5,   # Above max y for group A
            group == "B" ~ median_y,      # Centered on median y for group B
            group == "C" ~ median_y - 0.5 # Below median y for group C
        )
    )

# Create the base plot
savename <- "test.pdf"
pdf(savename, height=25, width =25)
plots <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_point() +
    # Add geom_text() with custom positioning by group
    geom_text(
        data = position_data,
        aes(x = text_x, y = text_y, label = label, color = group),
         size = 4, hjust = 0, vjust = 0
    ) +
    theme_minimal() +
    labs(titountle = "Custom Positioning of Text by Group in ggplot2")
plots
dev.off()
