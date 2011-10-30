setwd("/Workspace/R") # Working directory
library(foreign) # set up read DBF files

# SRF calc function using Jolly II for unequal sample units. 
# This needs a summary table, with each row a single transect:
## RSOD: number seen (COUNT).
## TRAN: area of sample (STRIP_AREA), area of transect (SUBUNIT_AR)
jolly.srf <- function(srftable) {
	obs <- sum(srftable$COUNT)
	R <- sum(srftable$COUNT)/sum(srftable$STRIP_AREA) # density
	Y <- R * sum(srftable$SUBUNIT_AR) # Pop estimate
	n <- length(srftable$Transect) # number of samples
	sy2 <- var(srftable$COUNT) #(1/(n-1))*(sum(rsod.sub$COUNT^2)-(sum(rsod.sub$COUNT)^2/n))            # Use var
	sz2 <- var(srftable$STRIP_AREA) #(1/(n-1))*(sum(tran$STRIP_AREA^2)-(sum(tran$STRIP_AREA)^2/n))
	szy <- (1/(n-1)) * (sum(srftable$COUNT*srftable$STRIP_AREA) - sum(srftable$SUBUNIT_AR) * sum(srftable$COUNT)/n)
	N <- length(srftable$Transect) * sum(srftable$SUBUNIT_AR) / sum(srftable$STRIP_AREA)	# Population of samples (NOT the real one - estimate!)
	varY <- (N*(N-n)/n)*sy2-2*R*n+R^2*sz2	# Variance of population estimate
	sterr <- sqrt(varY)
	answer <- c(obs, R, Y, n, N, varY, sterr)
	return(answer)
}

survey='surv01' # Define survey code.

# Get tran and rsod summary tables from csv:
fsod <- read.csv("/workspace/R/Test_Zam09data_FSOD.csv")
rsod <- read.csv("/workspace/R/Test_Zam09data_RSOD.csv")

strata <- levels(fsod$STRATUM) # list of all strata
spplist <- levels(rsod$CODE) #all species
adminlist <- levels(fsod$ADMIN) # all admin areas 

# housekeeping, set up results table.
framesize <- length(strata) * length(spplist) * length(adminlist) # maximum possible results table rowsize.
stratum <- strata; admin <- adminlist; sp <- spplist; obs <- numeric(framesize); R <- numeric(framesize); Y <- numeric(framesize); n <- numeric(framesize); N <- numeric(framesize); varY <- numeric(framesize); sterr <- numeric(framesize) # this sets the proper length of all vectors to make the results table.
srf.results <- data.frame(stratum, admin, sp, obs, R, Y, n, N, varY, sterr) # data frame with max required rows, correct levels and vector names.
result <- 1 # the row at which to start adding results to srf.results

# Here's the loop and calculation ...
system.time({  # To track how long it takes to do this loop.
for (stratum in strata) {
	#	fsod.substrat <- subset(fsod, fsod$STRATUM==stratum)                          
	fsod.substrat <- fsod[fsod$STRATUM == stratum,] 
	adminlist <- levels(factor(fsod.substrat$ADMIN)[, drop = TRUE])                         # SHould drop unused levels...
	#	rsod.substrat <- subset(rsod, rsod$STRATUM==stratum)
	rsod.substrat <- rsod[rsod$STRATUM == stratum,] 
	for (admin in adminlist) {
		#rsod.subadmin <- subset(rsod.substrat, rsod.substrat$ADMIN==admin)           # here too...
		rsod.subadmin <- rsod.substrat[rsod.substrat$ADMIN == admin,]
		spplist <- levels(factor(rsod.subadmin$CODE)[, drop = TRUE]) # rewrite this to remove COUNT=0
		# fsod.sub <- subset(fsod.substrat, fsod.substrat$ADMIN==admin) # subset of FSO data      AND HERE
		fsod.sub <- fsod.substrat[fsod.substrat$ADMIN==admin,]
		for (sp in spplist) {
			# rsod.subsp <- subset(rsod.subadmin, rsod.subadmin$CODE==sp) # make a subset of RSO data with only data for species sp 
			rsod.subsp <- rsod.subadmin[rsod.subadmin$CODE==sp, ]
			if (length(rsod.subsp$CODE) == 0) break
			rsod.tran <- aggregate(rsod.subsp[c(6)], by = list(Transect=rsod.subsp$TRAN), sum) # this gives only the transects that exist in rsod, but we need even the zeroes.
			if(length(fsod.sub$STRATUM) == 0 || is.na(fsod.sub$STRATUM[1])) break
			tran <- aggregate(fsod.sub[, c(18,19)], by=list(Transect=fsod.sub$TRAN), sum) # summary function to add together tranarea / samplearea, sorted by transects
			if (sum(tran$TRANAREA==0) ) break
			rsod.tran2 <- merge(tran, rsod.tran, all.x = TRUE) # gives us a table with area and count info, NAs for zero counts.
			rsod.tran2$COUNT[(is.na(rsod.tran2$COUNT))] <- 0 # clean up NA values
			srf.result <- jolly.srf(rsod.tran2) # ... and do the calc on species 
			srf.result.row <- srf.results[1,]
			srf.result.row[1] <- stratum
			srf.result.row[2] <- admin
			srf.result.row[3] <- sp
			srf.result.row[4:10] <- srf.result
			srf.results[result, ] <- srf.result.row # add the result
			result <- result + 1
		}
	}
}

}) # ends system timing.

srf.results <- srf.results[1:result-1, ] # trims the table to the actual number of results
write.csv(srf.results, paste("/workspace/Databases/test_", survey, "_result.csv", sep = ""))

# let the user know it's all done.
message <- 	paste("Finished. Table saved to: /Workspace/Databases/test_", survey, "_result.dbf", sep = "")
print(message)

