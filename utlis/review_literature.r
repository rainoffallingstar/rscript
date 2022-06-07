# libraries
library(rcoreoa)
library(RCurl)

# keys
CORE_KEY      <- "dBD1teHoYyjV6FsGQuJhC2wTAbm7L0lrP"
IBM_DISC_KEY  <- ""

# frames
QUERY    <- data.frame("all_of_the_words" = "algorithmic trading",
                       "find_those_words" = "b",       # fetch in titles
                       "language" = "en",
                       "year_from" = 2017,
                       "year_to" = 2020)

RES      <-  data.frame(ID            = character(),
                        Author_1st    = as.character(),
                        Title         = as.character(), 
                        Year          = as.character(), 
                        Description   = as.character(), 
                        Download_link = as.character())


# initial query counting how much pages of 100 articles exists
total_pages <- round( 
                      rcoreoa::core_advanced_search(
                                query = QUERY, 
                                key = CORE_KEY)
                      $totalHits / 100
                      )
total_pages <- 8

for (pp in 1:total_pages){
  res <- rcoreoa::core_advanced_search(query = QUERY, 
                                       key = CORE_KEY,
                                       page = pp,
                                       limit = 100)
  
  RES <- rbind(RES,
              data.frame(ID            = res$data$`_source`$id,
                         Author_1st    = unlist(lapply(res$data$`_source`$authors, function(x) x[1])),
                         Title         = res$data$`_source`$title, 
                         Year          = res$data$`_source`$year, 
                         Description   = unlist(res$data$`_source`$description), 
                         Download_link = res$data$`_source`$downloadUrl)
  )
}

# remove records without download link or description
who <- which(RES$Download_link != "" & RES$Description != "")
RES <- RES[who,]


# save auxiliar csv
write.csv2(RES,"c:/dev/algo_trading.csv", row.names = FALSE)
RES = read.csv2("c:/dev/algo_trading.csv")


download.file(
  url = as.vector(RES$Download_link),
  destfile = as.vector(get_pdf_name(RES)),
  mode = "wb"
)

for (i in 1:nrow(RES)){
  my_destfile <- paste("C:/dev/articles/",
                        RES$Year[i],
                        " ",
                        RES$ID[i],
                        " ", 
                        RES$Title[i],
                        ".pdf", sep="")
  
  my_url <-  as.character(RES$Download_link[i])
  print(c(my_url, my_destfile))
  
  tryCatch(download.file(url = my_url, 
                         destfile = my_destfile,
                         mode = "wb"), 
           error = function(e) e, 
           finally = print("Downloaded "))

}
