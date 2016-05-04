library(httr)
library(jsonlite)

CREEDS_URL <- 'http://amp.pharm.mssm.edu/CREEDS/'
response <- GET(paste0(CREEDS_URL, 'api'), query=list(id ='gene:84'))
if (response$status_code == 200){
	response <- fromJSON(httr::content(response, 'text'))
	print(response)
}
