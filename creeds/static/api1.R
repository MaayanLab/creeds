library(httr)
library(jsonlite)

CREEDS_URL <- 'http://amp.pharm.mssm.edu/CREEDS/'
response <- GET(paste0(CREEDS_URL, 'search'), query=list(q ='STAT3'))
if (response$status_code == 200){
	response <- fromJSON(httr::content(response, 'text'))
	print(response)
}
