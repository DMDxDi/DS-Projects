# Import libraries
pacman :: p_load(pacman,stats, dplyr, knitr, ggplot2, plotly)

# Read data
data <- read.table("/Users/ddx/R Scripts/EDA/EDA_Data Science Salaries in 2023/ds_salaries.csv", sep = ",", header = T)

# Nr. of instances
cat("Number of instances:", nrow(data))

# View data
head(data, 10)

       ### Data Cleaning ###

# Inspect data
check <- function(data) {
  l <- list()
  columns <- names(data)
  for (col in columns) {
    instances <- sum(!is.na(data[[col]]))
    dtypes <- class(data[[col]])
    unique <- length(unique(data[[col]]))
    sum_null <- sum(is.na(data[[col]]))
    duplicates <- sum(duplicated(data))
    l[[length(l) + 1]] <- c(col, dtypes, instances, unique, sum_null, duplicates)
  }
  data_check <- as.data.frame(do.call(rbind, l))
  names(data_check) <- c("column", "dtype", "instances", "unique", "sum_null", "duplicates")
  return(data_check)
}

check(data)

str(df)
describe(df)

                              ### EDA ###

      # Analysis for work_year #

wy_categ <- as.factor(ifelse(df$work_year  == 2020, '2020',
                             ifelse(df$work_year == 2021, '2021', 
                                    ifelse(df$work_year == 2022, '2022', 
                                           ifelse(df$work_year== 2023, '2023', '2020')))))

options(repr.plot.width=16, repr.plot.height=8)
my_palette <- c("#F8EDED", "#F6DFEB", "#E4BAD4", "#CE97B0")

wy_barchart <- ggplot(data.frame(wy_categ), aes(x = wy_categ)) +
  geom_bar(aes(fill = wy_categ))  +
  scale_fill_manual(values = my_palette) +
  ggtitle("Bar Chart for Work Year") +
  xlab("Year") +
  ylab("Frequency") +
  labs(fill = "Year") +
  stat_count(geom = "text", aes(label = after_stat(count)), vjust = -0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(color = "black", size = 20, face = "bold"),
    plot.subtitle = element_text(color = "#F6CD90",size = 12, face = "bold"),
    plot.caption = element_text(face = "italic"))

data <- data.frame(
  category=c("2020", "2021", "2022", "2023"),
  count=c(76, 230, 1664, 1785)
)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(round(data$fraction*100, 1), "%")

wy_piechart <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  ggtitle("Pie Chart for Work Year") +
  scale_fill_manual(values = my_palette)  +
  coord_polar(theta="y") +
  theme_void() +
  labs(fill = "Year") +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"))

grid.arrange(wy_barchart, wy_piechart, ncol = 2)

      # Analysis for experience_level #

# table(df$experience_level)

el_categ <- as.factor(ifelse(df$experience_level == 'EN', 'Entry-level',
                             ifelse(df$experience_level == 'MI', 'Mid-level', 
                                    ifelse(df$experience_level == 'SE', 'Senior-level', 
                                           ifelse(df$experience_level == 'EX', 'Executive-level', '')))))

options(repr.plot.width=16, repr.plot.height=8)
my_palette <- c("#FFF2F2", "#E5E0FF", "#8EA7E9", "#7286D3")

el_barchart <- ggplot(data.frame(el_categ), aes(x = el_categ)) +
  geom_bar(aes(fill = el_categ))  +
  scale_fill_manual(values = my_palette) +
  ggtitle("Bar Chart for Experience Level") +
  xlab("Level") +
  ylab("Frequency") +
  labs(fill = "Level") +
  stat_count(geom = "text", aes(label = after_stat(count)), vjust = -0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"),
    plot.subtitle = element_text(color = "#F6CD90",size = 12, face = "bold"),
    plot.caption = element_text(face = "italic"))

data <- data.frame(
  category=c("Entry-level", "Mid-level", "Senior-level", "Executive-level"),
  count=c(320, 805, 2516, 114)
)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(round(data$fraction*100, 1), "%")

el_piechart <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  ggtitle("Pie Chart for Experience Level") +
  scale_fill_manual(values = my_palette)  +
  coord_polar(theta="y") +
  theme_void() +
  labs(fill = "Level") +
  theme(
    plot.title = element_text(color = "black", size = 20, face = "bold"))

grid.arrange(el_barchart, el_piechart, ncol = 2)

      # Analysis for job_title #

# table(df$job_title)
cat("The dataset grasps", length(unique(df$job_title))," distinct job titles.")

# Compute the top 20 job titles in descending order
top20_job_titles <- head(sort(table(df$job_title), decreasing = T), 20)

# Create a bar plot with plotly
plot_ly(x = top20_job_titles, y = names(top20_job_titles), type = "bar",
        text = top20_job_titles, orientation = 'h', 
        marker = list(color = "#E1CCEC")) %>%
  layout(title = "Top 20 Jobs in Data Science", xaxis = list(title = "Count"), 
         yaxis = list(title = "Job Title")) 


      # Analysis for salary_in_usd #

options(scipen = 999)

# Central Tendencies
# Mean
cat("Mean:", mean(df$salary_in_usd))
cat("\n")

# Median
cat("Median:", median(df$salary_in_usd))
cat("\n")

# Mode
mode <- function(x){
  ta <- table(x)
  tam <- max(ta)
  if(all(ta==tam))
    mod <- NA
  else
    if(is.numeric(x))
      mod <- as.numeric(names(ta)[ta==tam])
  else
    mod <- names(ta)[ta==tam]
  return(mod)
}

cat("Mode:", mode(df$salary_in_usd))
cat("\n")

# Measure of Variability

# Std. Dev
cat("Standard Deviation:", sd(df$salary_in_usd))
cat("\n")

# Variance
cat("Variance:", var(df$salary_in_usd))
cat("\n")

# IQR
cat("Interquartile Range:") 
quantile(df$salary_in_usd)

options(repr.plot.width=18, repr.plot.height=6) 
require(gridExtra)

salary_hist <- ggplot(df, aes(x = salary_in_usd)) +
  geom_histogram(color = '#5c082c', fill ='#F5B0CB', bins=30) +
  labs(title = "Histogram for Salary ($)",x = "Salary",y = "Count") +
  theme_classic() +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"),
    plot.subtitle = element_text(color = "#F5B0CB",size = 12, face = "bold"),
    plot.caption = element_text(face = "italic"))

salary_boxplot <- ggplot(df, aes_string(x = df$salary_in_usd)) +
  geom_boxplot(outlier.colour = "#F5B0CB", outlier.shape = 11, outlier.size = 2, col = "#5c082c", notch = F) +
  labs(title = "Box Plot for Salary ($)",x = "Salary") +
  theme_classic() +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"))


grid.arrange(salary_hist,salary_boxplot, ncol=2)

plot(density(df$salary_in_usd),
     col="#F5B0CB",
     main="Density Plot for Salary",
     xlab="Salary",
     ylab="Density")
polygon(density(df$salary_in_usd),
        col="#F5B0CB")

      # Analysis for employment_type #


# table(df$employment_type)

et_categ <- as.factor(ifelse(df$employment_type == "CT", "Contract",
                             ifelse(df$employment_type == "FL", "Freelance", 
                                    ifelse(df$employment_type == "FT", "Full-Time", 
                                           ifelse(df$employment_type == "PT", "Part-Time", "")))))

options(repr.plot.width=16, repr.plot.height=8)
my_palette <- c("#FEDEFF", "#93C6E7", "#AEE2FF", "#B9F3FC")

et_barchart <- ggplot(data.frame(et_categ), aes(x = et_categ)) +
  geom_bar(aes(fill = et_categ))  +
  scale_fill_manual(values = my_palette) +
  ggtitle("Bar Chart for Employment Type") +
  xlab("Type") +
  ylab("Frequency") +
  labs(fill = "Level") +
  stat_count(geom = "text", aes(label = after_stat(count)), vjust = -0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"),
    plot.subtitle = element_text(color = "#F6CD90",size = 12, face = "bold"),
    plot.caption = element_text(face = "italic"))

data <- data.frame(
  category=c("Contract", "Freelance", "Full-Time", "Part-Time"),
  count=c(10, 10, 3718, 17)
)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(round(data$fraction*100, 1), "%")

et_piechart <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  ggtitle("Pie Chart for Employment Type") +
  scale_fill_manual(values = my_palette)  +
  coord_polar(theta="y") +
  theme_void() +
  labs(fill = "Type") +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"))

grid.arrange(et_barchart, et_piechart, ncol = 2)

# Analysis for remote_ratio #
# table(df$remote_ratio)

rr_categ <- as.factor(ifelse(df$remote_ratio  == 0, "On-site",
                             ifelse(df$remote_ratio == 50, "Hybrid", 
                                    ifelse(df$remote_ratio == 100, "Remote", ""))))

options(repr.plot.width=16, repr.plot.height=8)
my_palette <- c("#f0dfd1", "#f0c9a8", "#edb482")

rr_barchart <- ggplot(data.frame(rr_categ ), aes(x = rr_categ)) +
  geom_bar(aes(fill = rr_categ ))  +
  scale_fill_manual(values = my_palette) +
  labs(fill = "Type") +
  stat_count(geom = "text", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Bar Chart for  % Remote Work",x = "Type", y="Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"),
    plot.subtitle = element_text(color = "#F6CD90",size = 12, face = "bold"),
    plot.caption = element_text(face = "italic"))

data <- data.frame(
  category=c("On-site", "Hybrid", "Remote"),
  count=c(1923, 189, 1643)
)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(round(data$fraction*100, 1), "%")

rr_piechart <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  ggtitle("Pie Chart for % Remote Work") +
  scale_fill_manual(values = my_palette)  +
  coord_polar(theta="y") +
  theme_void() +
  labs(fill = "Type") +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"))

grid.arrange(rr_barchart, rr_piechart, ncol = 2)

      # Analysis for company_size #
# table(df$company_size)

cs_categ <- as.factor(ifelse(df$company_size  == "L", "Large",
                             ifelse(df$company_size == "M", "Medium", 
                                    ifelse(df$company_size == "S", 'Small', ''))))

options(repr.plot.width=16, repr.plot.height=8)
my_palette <- c("#C8E3D4", "#96C7C1", "#89B5AF")

cs_barchart <- ggplot(data.frame(cs_categ ), aes(x = cs_categ)) +
  geom_bar(aes(fill = cs_categ ))  +
  scale_fill_manual(values = my_palette) +
  ggtitle("Bar Chart for Companies Size") +
  xlab("Size") +
  ylab("Frequency") +
  labs(fill = "Type") +
  stat_count(geom = "text", aes(label = after_stat(count)), vjust = -0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"),
    plot.subtitle = element_text(color = "#F6CD90",size = 12, face = "bold"),
    plot.caption = element_text(face = "italic"))

data <- data.frame(
  category=c("Large", "Medium", "Small"),
  count=c(454, 3153, 148)
)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(round(data$fraction*100, 1), "%")

cs_piechart <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  ggtitle("Pie Chart for Companies Size") +
  scale_fill_manual(values = my_palette)  +
  coord_polar(theta="y") +
  theme_void() +
  labs(fill = "Size") +
  theme(
    plot.title = element_text(color = "#383335", size = 20, face = "bold"))

grid.arrange(cs_barchart, cs_piechart, ncol = 2)


      # Analysis for company_location #

# table(df$company_location)

cat("The dataset grasps", length(unique(df$company_location))," distinct company locations.")

# Compute the top 20 company locations in descending order
top20_company_location <- head(sort(table(df$company_location), decreasing = T), 20)

# Create a bar plot with plotly
plot_ly(x = top20_company_location, y = names(top20_company_location), type = "bar",
        text = top20_job_titles, orientation = 'h', 
        marker = list(color = "#DDB6C6")) %>%
  layout(title = "Top 20 Data Science Company Locations", xaxis = list(title = "Count"), 
         yaxis = list(title = "Company Location")) 


      # Analysis for employee_residence #

# table(df$employee_residence)

cat("The dataset grasps", length(unique(df$employee_residence))," distinct employees residence.")

# Compute the top 20 employees residence in descending order
top20_employee_residence <- head(sort(table(df$employee_residence), decreasing = T), 20)

# Create a bar plot with plotly
plot_ly(x = top20_employee_residence, y = names(top20_employee_residence), type = "bar",
        text = top20_job_titles, orientation = 'h', 
        marker = list(color = "#A6E3E9")) %>%
  layout(title = "Top 20 Data Science Employees Residence", xaxis = list(title = "Count"), 
         yaxis = list(title = "Employee Residence"))


      # Analysis for Salary by Year #
options(scipen = 999)

options(repr.plot.width=20, repr.plot.height=10) 
require(gridExtra)

my_palette <- c("#F8EDED", "#F6DFEB", "#E4BAD4", "#CE97B0")

wy_sl_boxplot <- ggplot(df, aes(x=wy_categ, y=salary_in_usd)) + 
  geom_boxplot(fill=my_palette,
               outlier.colour = "#FFB9B9", 
               outlier.shape = 11, outlier.size = 2, 
               col = my_palette, 
               notch = F) +
  labs(
    title = "Distribution of Salaries by Year",
    x= "Year",
    y="Salary($)") 


wy_sl_hist <- ggplot(df, aes(x=salary_in_usd, fill=wy_categ)) + 
  geom_histogram(color="black", bins = 30) +
  labs(
    title = "Distribution of Salaries by Year",
    x="Salary", 
    y="Frequency",
    fill="Year") +
  scale_fill_manual(values = c("#F8EDED", "#F6DFEB", "#E4BAD4", "#CE97B0")) 


grid.arrange(wy_sl_boxplot, wy_sl_hist, ncol=2)


      # Analysis for Salary by Experience Level #
options(scipen = 999)

options(repr.plot.width=20, repr.plot.height=10) 
require(gridExtra)

my_palette <- c("#FFF2F2", "#E5E0FF", "#8EA7E9", "#7286D3")

el_sl_boxplot <- ggplot(df, aes(x=el_categ, y=salary_in_usd)) + 
  geom_boxplot(fill=my_palette,
               outlier.colour = "#B2A4FF", 
               outlier.shape = 11, outlier.size = 2, 
               col = my_palette, 
               notch = F) +
  labs(
    title = "Distribution of Salaries by Experience Level",
    x= "Level",
    y="Salary($)") 


el_sl_hist <- ggplot(df, aes(x=salary_in_usd, fill=el_categ)) + 
  geom_histogram(color="black", bins = 30) +
  labs(
    title = "Distribution of Salaries by Experience Level",
    x="Salary", 
    y="Frequency",
    fill="Level") +
  scale_fill_manual(values = c("#FFF2F2", "#E5E0FF", "#8EA7E9", "#7286D3")) 


grid.arrange(el_sl_boxplot, el_sl_hist, ncol=2)


      # Analysis for Salary by Employment Type #
options(scipen = 999)

options(repr.plot.width=20, repr.plot.height=10) 
require(gridExtra)

my_palette <- c("#FEDEFF", "#93C6E7", "#AEE2FF", "#B9F3FC")

et_sl_boxplot <- ggplot(df, aes(x=et_categ, y=salary_in_usd)) + 
  geom_boxplot(fill=my_palette,
               outlier.colour = "#FFF2F2", 
               outlier.shape = 11, outlier.size = 2, 
               col = my_palette, 
               notch = F) +
  labs(
    title = "Distribution of Salaries by Employment Type",
    x= "Type",
    y="Salary($)")

et_sl_hist <- ggplot(df, aes(x=salary_in_usd, fill=et_categ)) + 
  geom_histogram(color="black", bins = 30) +
  labs(
    title = "Distribution of Salaries by Employment Type",
    x="Salary", 
    y="Frequency",
    fill="Type") +
  scale_fill_manual(values = c("#FEDEFF", "#93C6E7", "#AEE2FF", "#B9F3FC")) 

grid.arrange(et_sl_boxplot, et_sl_hist, ncol=2)


      # Analysis for Salary by Company Size #
options(scipen = 999)

options(repr.plot.width=20, repr.plot.height=10) 
require(gridExtra)

my_palette <- c("#C8E3D4", "#96C7C1", "#89B5AF")

cs_sl_boxplot <- ggplot(df, aes(x=cs_categ, y=salary_in_usd)) + 
  geom_boxplot(fill=my_palette,
               outlier.colour = "#305F72", 
               outlier.shape = 11, outlier.size = 2, 
               col = my_palette, 
               notch = F) +
  labs(
    title = "Distribution of Salaries by Company Size",
    x= "Size",
    y="Salary($)") 

cs_sl_hist <- ggplot(df, aes(x=salary_in_usd, fill=cs_categ)) + 
  geom_histogram(color="black", bins = 30) +
  labs(
    title = "Distribution of Salaries by Company Size",
    x="Salary", 
    y="Frequency",
    fill="Size") +
  scale_fill_manual(values = c("#C8E3D4", "#96C7C1", "#89B5AF")) 

grid.arrange(cs_sl_boxplot, cs_sl_hist, ncol=2)


      # Analysis for Salary by % Remote Work #
options(scipen = 999)

options(repr.plot.width=20, repr.plot.height=10) 
require(gridExtra)

my_palette <- c("#f0dfd1", "#f0c9a8", "#edb482")

rr_sl_boxplot <- ggplot(df, aes(x=rr_categ, y=salary_in_usd)) + 
  geom_boxplot(fill=my_palette,
               outlier.colour = "#FFD966", 
               outlier.shape = 11, outlier.size = 2, 
               col = my_palette, 
               notch = F) +
  labs(
    title = "Distribution of Salaries by % Remote Work",
    x= "Type",
    y="Salary($)")

rr_sl_hist <- ggplot(df, aes(x=salary_in_usd, fill=rr_categ)) + 
  geom_histogram(color="black", bins = 30) +
  labs(
    title = "Distribution of Salaries by % Remote Work",
    x="Salary", 
    y="Frequency",
    fill="Type") +
  scale_fill_manual(values = c("#f0dfd1", "#f0c9a8", "#edb482")) 

grid.arrange(rr_sl_boxplot, rr_sl_hist, ncol=2)






