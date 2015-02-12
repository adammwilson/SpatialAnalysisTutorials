
          
          library(getopt)
          ## get options
          opta <- getopt(matrix(c(
          'date', 'd', 1, 'character'
          ), ncol=4, byrow=TRUE))
          ## extract value
          date=as.Date(opta$date) 
          
          ## Now your script using date as an input
          
          print(date+1)
          q("no")
          