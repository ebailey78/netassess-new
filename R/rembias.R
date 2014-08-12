rembias <- function() 
{
    options(warn = -1)
    require(sp)
    require(deldir)
    require(rgdal)
    require(maptools)
    require(fgui)
    require(tcltk)
    require(sqldf)
    pb <- winProgressBar(title = "Removal Bias Tool Starting Up", 
        min = 0, max = 1, width = 500, label = "Loading data. Please wait")
    for (i in j <- c("ozone", "pm25avg", "contpm25_all", "latlongs", 
        "o3dvs", "pm25dvs", "contpm25dvs", "addresses")) {
        data(list = i, package = "netassess")
        setWinProgressBar(pb, grep(i, j)/length(j))
    }
    Sys.sleep(2)
    close(pb)
    tolambert <- function(xycoords) {
        tempxy <- SpatialPoints(xycoords, proj4string = CRS("+proj=longlat +ellps=clrk66"))
        tempxy <- spTransform(tempxy, CRS("+proj=lcc +lon_0=97w +lat_0=40n +lat_1=33n +lat_2=45n +ellps=clrk66"))
        dimnames(tempxy@coords)[[2]] <- c("lamx", "lamy")
        return(tempxy@coords)
    }
    choose.poll <- function() {
        if (grepl("ozone", varin)) {
            crit.val <<- "maxval"
            dvs <<- o3dvs
            poll <- ozone[as.numeric(substr(ozone$colldate, 6, 
                7)) %in% 5:9 & ozone$yr == 2008, ]
            compsites <- sqldf("select distinct siteid\n\t\t\t\t\t\t\t\tfrom (select distinct siteid, count(*) as ndays \n\t\t\t\t\t\t\t\t\t\t\tfrom poll where (obspct>=75 or maxval>0.075)\n\t\t\t\t\t\t\t\t\t\t\tgroup by siteid) \n\t\t\t\t\t\t\t\twhere ndays>=153*0.75", 
                method = "raw")$siteid
        }
        else if (grepl("pm25frm|contpm25", varin)) {
            if (grepl("contpm25", varin)) {
                crit.val <<- "contpm25"
                dvs <<- contpm25dvs
                poll <- contpm25_all[as.numeric(contpm25_all$yr) == 
                  2008, ]
                daysuse <- as.character(seq(from = as.Date("2005-01-01", 
                  "%Y-%m-%d"), to = as.Date("2008-12-31", "%Y-%m-%d"), 
                  by = 1))
            }
            else {
                crit.val <<- "pm25"
                dvs <<- pm25dvs
                poll <- pm25avg[as.numeric(pm25avg$yr) == 2008, 
                  ]
                if (grepl("3-day", varin)) {
                  sampsched <- "3-day"
                }
                else {
                  sampsched <- "6-day"
                }
                startday <- list("2005-01-01", "2005-01-04")
                names(startday) <- c("3-day", "6-day")
                daysuse <- as.character(seq(from = as.Date(startday[[sampsched]], 
                  "%Y-%m-%d"), to = as.Date("2008-12-31", "%Y-%m-%d"), 
                  by = as.numeric(substr(sampsched, 1, 1))))
            }
            poll$qtr <- ceiling(as.numeric(substr(poll$colldate, 
                6, 7))/3)
            daysuse <- daysuse[as.numeric(substr(daysuse, 1, 
                4)) == 2008]
            daysuse <- data.frame(days = daysuse, yr = as.numeric(substr(daysuse, 
                1, 4)), qtr = ceiling(as.numeric(substr(daysuse, 
                6, 7))/3), stringsAsFactors = FALSE)
            ndays <- sqldf("select distinct qtr, count(*) as ndays\n\t\t\t\tfrom daysuse\n\t\t\t\tgroup by qtr", 
                method = "raw")
            compsites <- sqldf("select distinct siteid, qtr, count(*) as nobs\n\t\t\t\tfrom poll\n\t\t\t\tgroup by siteid, qtr", 
                method = "raw")
            compsites <- merge(compsites, ndays, by = "qtr", 
                all.x = TRUE)
            compsites <- sqldf("select distinct siteid\n\t\t\t\t\t\t\t\t\t\t\tfrom (select distinct siteid, count(*) as nqtrs\n\t\t\t\t\t\t\t\t\t\t\t\t\t\tfrom compsites\n\t\t\t\t\t\t\t\t\t\t\t\t\t\twhere nobs>=0.75*ndays\n\t\t\t\t\t\t\t\t\t\t\t\t\t\tgroup by siteid)\n\t\t\t\t\t\t\t\t\t\t\twhere nqtrs=4", 
                method = "raw")$siteid
        }
        poll <<- poll[poll$siteid %in% compsites, ]
    }
    calc.rembias <- function() {
        done <<- 3
        if (tclvalue(f.but$object) != "" & length(test.sites) == 
            0) {
            test.sites <- readLines(tclvalue(f.but$object))
            test.sites <<- unlist(strsplit(test.sites, ",", fixed = TRUE))
            choose.poll()
        }
        byvars <- split(poll, poll$colldate)
        pb <- winProgressBar(title = "Calculating Bias", label = "0% done", 
            max = 100)
        count <<- 0
        rembias <- lapply(byvars, function(x) {
            garbage2 <- lapply(test.sites, function(z) {
                count <<- count + 1
                x <- merge(x, latlongs, by = "siteid", all.x = TRUE)
                x <- cbind(x, tolambert(x[, c("longitude", "latitude")]))
                useset <- x[!x$siteid %in% test.sites, ]
                testset <- x[x$siteid %in% test.sites, ]
                delsgs <- data.frame(deldir(x = useset$longitude, 
                  y = useset$latitude, dpl = list(x = testset$longitude[testset$siteid == 
                    z], y = testset$latitude[testset$siteid == 
                    z]))$delsgs)
                delsgs <- delsgs[delsgs$ind1 == max(delsgs$ind1) | 
                  delsgs$ind2 == max(delsgs$ind1), ]
                sitesused <- unique(c(delsgs$ind1, delsgs$ind2))
                sitesused <- sitesused[sitesused != max(sitesused)]
                sitesused <- list(useset[sitesused, "siteid"])
                delsgs$distance <- with(delsgs, sqrt((x1 - x2)^2 + 
                  (y1 - y2)^2))
                delsgs$weights <- with(delsgs, (1/distance^2)/sum(1/distance^2))
                delsgs <- data.frame(ind2 = with(delsgs, ifelse(ind1 == 
                  max(ind1), ind2, ind1)), weights = delsgs$weights, 
                  stringsAsFactors = FALSE)
                delsgs <- data.frame(conc = useset[delsgs$ind2, 
                  crit.val], weights = delsgs$weights, stringsAsFactors = FALSE)
                testset$intconc[testset$siteid == z] <- round(with(delsgs, 
                  sum(conc * weights)), digits = 3)
                pbval <- round(count/(length(byvars) * length(test.sites)) * 
                  100, digits = 0)
                setWinProgressBar(pb, value = pbval, title = "Calculating Bias", 
                  label = paste(pbval, "% done", sep = ""))
                testset$sitesused <- sitesused
                return(testset[testset$siteid == z, ])
            })
            garbage2 <- do.call("rbind", garbage2)
            return(garbage2)
        })
        Sys.sleep(0.5)
        close(pb)
        rembias <- do.call("rbind", rembias)
        rembias$bias <- with(rembias, eval(parse(text = paste("intconc-", 
            crit.val, sep = ""))))
        rembias <- by(rembias, list(rembias$siteid), function(x) {
            temp <- data.frame(siteid = unique(x$siteid), meanbias = mean(x$bias, 
                na.rm = TRUE), sdbias = sd(x$bias, na.rm = TRUE), 
                nobs = NROW(x$bias), stringsAsFactors = FALSE)
            sitetable <- table(unlist(x$sitesused))
            temp$sitesused <- list(paste(names(sitetable), " (", 
                sitetable, ")", sep = ""))
            return(temp)
        })
        rembias <- do.call("rbind", rembias)
        rembias$pvalue <- pt(rembias$meanbias/rembias$sdbias * 
            sqrt(rembias$nobs), rembias$nobs - 1)
        rembias$pvalue <- ifelse(rembias$meanbias < 0, rembias$pvalue * 
            2, (1 - rembias$pvalue) * 2)
        rembias <- merge(rembias, latlongs, by = "siteid", all.x = TRUE)
        rembias <- cbind(rembias, tolambert(rembias[, c("longitude", 
            "latitude")]))
        if (grepl("ozone", varin)) {
            abs.breaks <- c(1e-10, 0.001, 0.003, 0.005, 0.01, 
                1e+10)
            rdig <- 3
        }
        else {
            abs.breaks <- c(1e-06, 0.5, 1, 1.5, 2, 1e+06)
            rdig <- 1
        }
        all.quantile <- c(-rev(abs.breaks), abs.breaks)
        col.palette <- colorRampPalette(c("blue", "white", "red"))(11)
        rembias$all.colors <- as.character(cut(rembias$meanbias, 
            breaks = unique(all.quantile), include.lowest = TRUE, 
            labels = col.palette))
        rembias$sig <- ifelse(rembias$pvalue > 0.05, "insig", 
            "sig")
        rembias$meanbias <- round(rembias$meanbias, digits = rdig)
        rembias <- merge(rembias, dvs, by = "siteid", all.x = TRUE)
        rembias <- merge(rembias, addresses, by = "siteid", all.x = TRUE)
        rembias <- rembias[, -grep("0507|lamx|lamy|address|nickname", 
            names(rembias))]
        remain.sites <- latlongs[latlongs$siteid %in% unique(poll$siteid[!poll$siteid %in% 
            test.sites]), ]
        remain.sites <- merge(remain.sites, addresses[, c("siteid", 
            "nickaddr")], by = "siteid", all.x = TRUE)
        remain.sites <- merge(remain.sites, dvs[, grep("siteid|0608", 
            names(dvs))], by = "siteid", all.x = TRUE)
        rembias <- merge(rembias, remain.sites, all = TRUE, sort = FALSE)
        coordpts <- SpatialPoints(rembias[, c("longitude", "latitude")], 
            proj4string = CRS("+proj=longlat +ellps=clrk66"))
        toshp <- SpatialPointsDataFrame(coordpts, rembias[, -grep("longitude|latitude|all.colors", 
            names(rembias))])
        toshp$sitesused <- sapply(toshp$sitesused, function(x) paste(x, 
            collapse = "  "))
        writePointsShape(toshp, paste(dname, "/", tclvalue(fname$object), 
            sep = ""))
        write.csv(data.frame(rembias[, -grep("all.colors|sitesused", 
            names(rembias))], sitesused = sapply(rembias$sitesused, 
            function(x) paste(x, collapse = " ")), stringsAsFactors = FALSE), 
            paste(dname, "/", tclvalue(fname$object), ".csv", 
                sep = ""), row.names = FALSE)
        kml <- file(paste(dname, "/", tclvalue(fname$object), 
            ".kml", sep = ""), open = "w")
        scale_file <- "rembias_pm25_scale.png"
        if (grepl("ozone", varin)) {
            parm <- "Ozone"
            scale_file <- "rembias_o3_scale.png"
        }
        else if (grepl("pm25frm", varin)) {
            parm <- "PM2.5 (FRM)"
        }
        else {
            parm <- "Continuous PM2.5"
        }
        cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", "<kml xmlns=\"http://earth.google.com/kml/2.1\">", 
            "<Document>", "\t<Folder>", paste("\t<name>Removed ", 
                parm, " Sites</name>", sep = ""), "\t<visibility>0</visibility>", 
            "\t<open>0</open>", "\t<ScreenOverlay>", "\t\t<overlayXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>", 
            "\t\t<screenXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>", 
            "\t\t<Icon>", paste("\t\t    <href>", scale_file, 
                "</href>", sep = ""), "\t\t\t<scale>1</scale>", 
            "\t\t</Icon>", "\t\t<ListStyle>", "\t\t\t<listItemType>checkHideChildren</listItemType>", 
            "\t\t</ListStyle>", "\t</ScreenOverlay>", sep = "\n", 
            file = kml)
        sigs <- split(rembias[rembias$sig == "sig", ], rembias$siteid[rembias$sig == 
            "sig"])
        insigs <- split(rembias[rembias$sig == "insig", ], rembias$siteid[rembias$sig == 
            "insig"])
        garbage <- lapply(sigs, function(x) {
            cat(paste("\t<Placemark>\n\n\t\t\t\t\t\t<name>", 
                x$siteid, "</name>\n\n\t\t\t\t\t\t<Style id=\"sigsites\">\n\n\t\t\t\t\t\t\t<IconStyle>\n\n\t\t\t\t\t\t\t\t<scale>1</scale>\n\n\t\t\t\t\t\t\t\t<color>ff", 
                substr(x$all.colors, 6, 7), substr(x$all.colors, 
                  4, 5), substr(x$all.colors, 2, 3), "</color>\n\n\t\t\t\t\t\t\t\t<Icon>\n\n\t\t\t\t\t\t\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/donut.png</href>\n\n\t\t\t\t\t\t\t\t</Icon>\n\n\t\t\t\t\t\t\t</IconStyle>\n\n\t\t\t\t\t\t</Style>\n\n\t\t\t\t\t\t<description>\n", 
                sep = ""), file = kml)
            if (grepl("ozone", varin)) {
                cat(paste("<![CDATA[<table width=\"300\"><tr><td><p>AQS ID: ", 
                  x$siteid, "<br>", x$nickaddr, "<br><br>2006-2008 Design Value: ", 
                  x$dv0608, " ppm (", round(x$dv0608/0.075 * 
                    100, digits = 0), "%)<br><br>Average Bias: ", 
                  x$meanbias, " ppm<br><br>Sites Used in Bias Calculation:<br>", 
                  sapply(x$sitesused, function(y) paste(y, collapse = "<br>")), 
                  "<br>(Numbers in parentheses reflect number of times site was used in bias calculation)</p></td></tr></table>]]>\n", 
                  sep = ""), file = kml)
            }
            else {
                cat(paste("<![CDATA[<table width=\"300\"><tr><td><p>AQS ID: ", 
                  x$siteid, "<br>", x$nickaddr, "<br><br>2006-2008 Annual Design Value: ", 
                  x$dva0608, " ug/m3 (", round(x$dva0608/15 * 
                    100, digits = 0), "%)<br>2006-2008 Daily Design Value: ", 
                  x$dvd0608, " ug/m3 (", round(x$dvd0608/35 * 
                    100, digits = 0), "%)<br><br>Average Bias: ", 
                  x$meanbias, " ug/m3<br><br>Sites Used in Bias Calculation:<br>", 
                  sapply(x$sitesused, function(y) paste(y, collapse = "<br>")), 
                  "<br>(Numbers in parentheses reflect number of times site was used in bias calculation)</p></td></tr></table>]]>\n", 
                  sep = ""), file = kml)
            }
            cat(paste("</description>\n\n\t\t\t\t\t\t<Point>\n\n\t\t\t\t\t\t\t<extrude>0</extrude>\n\n\t\t\t\t\t\t\t<altitudeMode>ClampToGround</altitudeMode>\n\n\t\t\t\t\t\t\t<coordinates>", 
                x$longitude, ",", x$latitude, ",0</coordinates>\n\n\t\t\t\t\t\t</Point>\n\n\t\t\t\t\t</Placemark>\n", 
                sep = ""), file = kml)
        })
        garbage <- lapply(insigs, function(x) {
            cat(paste("\t<Placemark>\n\n\t\t\t\t\t\t<name>", 
                x$siteid, "</name>\n\n\t\t\t\t\t<Style id=\"sigsites\">\n\n\t\t\t\t\t\t<IconStyle>\n\n\t\t\t\t\t\t\t<scale>1.5</scale>\n\n\t\t\t\t\t\t\t<color>ff", 
                substr(x$all.colors, 6, 7), substr(x$all.colors, 
                  4, 5), substr(x$all.colors, 2, 3), "</color>\n\n\t\t\t\t\t\t\t<Icon>\n\n\t\t\t\t\t\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n\n\t\t\t\t\t\t\t</Icon>\n\n\t\t\t\t\t\t</IconStyle>\n\n\t\t\t\t\t</Style>\n\n\t\t\t\t\t<description>\n", 
                sep = ""), file = kml)
            if (grepl("ozone", varin)) {
                cat(paste("<![CDATA[<table width=\"300\"><tr><td><p>AQS ID: ", 
                  x$siteid, "<br>", x$nickaddr, "<br><br>2006-2008 Design Value: ", 
                  x$dv0608, " ppm (", round(x$dv0608/0.075 * 
                    100, digits = 0), "%)<br><br>Average Bias: ", 
                  x$meanbias, " ppm<br><br>Sites Used in Bias Calculation:<br>", 
                  sapply(x$sitesused, function(y) paste(y, collapse = "<br>")), 
                  "<br>(Numbers in parentheses reflect number of times site was used in bias calculation)</p></td></tr></table>]]>\n", 
                  sep = ""), file = kml)
            }
            else {
                cat(paste("<![CDATA[<table width=\"300\"><tr><td><p>AQS ID: ", 
                  x$siteid, "<br>", x$nickaddr, "<br><br>2006-2008 Annual Design Value: ", 
                  x$dva0608, " ug/m3 (", round(x$dva0608/15 * 
                    100, digits = 0), "%)<br>2006-2008 Daily Design Value: ", 
                  x$dvd0608, " ug/m3 (", round(x$dvd0608/35 * 
                    100, digits = 0), "%)<br><br>Average Bias: ", 
                  x$meanbias, " ug/m3<br><br>Sites Used in Bias Calculation:<br>", 
                  sapply(x$sitesused, function(y) paste(y, collapse = "<br>")), 
                  "<br>(Numbers in parentheses reflect number of times site was used in bias calculation)</p></td></tr></table>]]>\n", 
                  sep = ""), file = kml)
            }
            cat(paste("</description>\n\n\t\t\t\t\t\t<Point>\n\n\t\t\t\t\t\t\t<extrude>0</extrude>\n\n\t\t\t\t\t\t\t<altitudeMode>ClampToGround</altitudeMode>\n\n\t\t\t\t\t\t\t<coordinates>", 
                x$longitude, ",", x$latitude, ",0</coordinates>\n\n\t\t\t\t\t\t</Point>\n\n\t\t\t\t\t</Placemark>\n", 
                sep = ""), file = kml)
        })
        cat("\t</Folder>\n\n\t\t\t<Folder>\n\n\t\t\t\t<Style id=\"noitems\">\n\n\t\t\t\t\t<ListStyle>\n\n\t\t\t\t\t\t<listItemType>checkHideChildren</listItemType>\n\n\t\t\t\t\t</ListStyle>\n\n\t\t\t\t</Style>\n\n\t\t\t\t<styleUrl>#noitems</styleUrl>\n\n\t\t\t\t<name>Remaining Sites</name>\n", 
            file = kml)
        byvars <- split(remain.sites, remain.sites$siteid)
        garbage <- lapply(byvars, function(x) {
            cat(paste("\t<Placemark>\n\n\t\t\t\t\t\t<name>", 
                x$siteid, "</name>\n\n\t\t\t\t\t<Style id=\"othrsites\">\n\n\t\t\t\t\t\t<IconStyle>\n\n\t\t\t\t\t\t\t<scale>0.75</scale>\n\n\t\t\t\t\t\t\t<color>ff808080</color>\n\n\t\t\t\t\t\t\t<Icon>\n\n\t\t\t\t\t\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n\n\t\t\t\t\t\t\t</Icon>\n\n\t\t\t\t\t\t</IconStyle>\n\n\t\t\t\t\t</Style>\n\n\t\t\t\t\t<description>\n", 
                sep = ""), file = kml)
            if (grepl("ozone", varin)) {
                cat(paste("<![CDATA[<table width=\"300\"><tr><td><p>AQS ID: ", 
                  x$siteid, "<br>", x$nickaddr, "<br><br>2006-2008 Design Value: ", 
                  x$dv0608, " ppm (", round(x$dv0608/0.075 * 
                    100, digits = 0), "%)</p></td></tr></table>]]>\n", 
                  sep = ""), file = kml)
            }
            else {
                cat(paste("<![CDATA[<table width=\"300\"><tr><td><p>AQS ID: ", 
                  x$siteid, "<br>", x$nickaddr, "<br><br>2006-2008 Annual Design Value: ", 
                  x$dva0608, " ug/m3 (", round(x$dva0608/15 * 
                    100, digits = 0), "%)<br>2006-2008 Daily Design Value: ", 
                  x$dvd0608, " ug/m3 (", round(x$dvd0608/35 * 
                    100, digits = 0), "%)</p></td></tr></table>]]>\n", 
                  sep = ""), file = kml)
            }
            cat(paste("</description>\n\n\t\t\t\t\t\t<Point>\n\n\t\t\t\t\t\t\t<extrude>0</extrude>\n\n\t\t\t\t\t\t\t<altitudeMode>ClampToGround</altitudeMode>\n\n\t\t\t\t\t\t\t<coordinates>", 
                x$longitude, ",", x$latitude, ",0</coordinates>\n\n\t\t\t\t\t\t</Point>\n\n\t\t\t\t\t</Placemark>\n", 
                sep = ""), file = kml)
        })
        cat("\t</Folder>\n\t\t\n\t\t</Document>\n\n\t\t</kml>", 
            sep = "\n", file = kml)
        close(con = kml)
        shell(paste("zip -j \"", gsub("/", "\\", dname, fixed = TRUE), 
            "\\", tclvalue(fname$object), ".kmz\" \"", gsub("/", 
                "\\", dname, fixed = TRUE), "\\", tclvalue(fname$object), 
            ".kml\" \"", gsub("/", "\\", .Library, fixed = TRUE), 
            "\\netassess\\", scale_file, "\"", sep = ""))
        shell(paste("del /q \"", gsub("/", "\\", dname, fixed = TRUE), 
            "\\", tclvalue(fname$object), ".kml\"", sep = ""))
        done <<- 2
    }
    data.test <- function(checkvar) {
        done <<- 0
        if (checkvar == "main3") {
            if ((length(test.sites) == 0 & tclvalue(f.but$object) == 
                "") | tclvalue(fname$object) == "" | dname == 
                "" | varin == "") {
                tt <- tkmessageBox(title = "'Doh!'", message = "One or more of the items marked by '*' is missing!", 
                  icon = "info", type = "ok")
                tcl("focus", "-force", dialog)
            }
            else {
                done <<- 1
            }
        }
    }
    newguiFilename <- function(sframe, text = "", default = "", 
        title = "", filter = "") {
        frame <- tkframe(sframe)
        var <- tclVar(default)
        filename <- tkentry(frame, width = 50, textvariable = var)
        cmd <- function() {
            tempstr <- tclvalue(tkgetOpenFile(title = title, 
                filetypes = filter))
            if (!nchar(tempstr)) 
                return()
            tclvalue(var) <- tempstr
        }
        file.but <- tkbutton(frame, text = text, command = cmd)
        tkgrid(file.but, filename)
        tkgrid.configure(file.but, sticky = "ew")
        return(list(object = var, guiObject = frame))
    }
    done <<- 0
    tclRequire("BWidget")
    while (done != 3) {
        dname <- ""
        test.sites <<- ""
        dialog <- tktoplevel()
        tkwm.title(dialog, "Removal Bias Tool")
        tkgrid(tklabel(dialog, text = "Valid choices in the tree below are colored \"blue\".\nYou may have to expand an individual branch\nof the tree by clicking on the \"plus\" to\nget to a valid choice. *", 
            justify = "left"), sticky = "nw")
        dframe <- guiFrame(sframe = dialog)
        tkgrid(dframe, sticky = "nw")
        scrx <- tkscrollbar(dframe, command = function(...) tkxview(polltree, 
            ...), orient = "horizontal")
        scry <- tkscrollbar(dframe, command = function(...) tkyview(polltree, 
            ...), orient = "vertical")
        polltree <- tkwidget(dframe, "Tree", xscrollcommand = function(...) tkset(scrx, 
            ...), yscrollcommand = function(...) tkset(scry, 
            ...))
        tkgrid(polltree, scry)
        tkgrid.configure(scry, sticky = "news")
        tkgrid(scrx, sticky = "ew")
        tkinsert(polltree, "end", "root", "ozone", text = "Ozone", 
            fill = "blue")
        tkinsert(polltree, "end", "root", "pm25frm", text = "PM2.5 (FRM)", 
            selectable = "no")
        tkinsert(polltree, "end", "pm25frm", "pm25frm3-day", 
            text = "3-day schedule", fill = "blue")
        tkinsert(polltree, "end", "pm25frm", "pm25frm6-day", 
            text = "6-day schedule", fill = "blue")
        tkinsert(polltree, "end", "root", "contpm25", text = "Continuous PM2.5", 
            fill = "blue")
        tcl(polltree, "bindText", "<ButtonRelease-1>", function(...) if (tclvalue(tcl(polltree, 
            "selection", "get")) != "") {
            varin <<- tclvalue(tcl(polltree, "selection", "get"))
        })
        dframe <- guiFrame(borderwidth = 0, sframe = dialog, 
            , relief = "flat")
        tkgrid.configure(dframe, sticky = "news")
        s.but <- tkbutton(dframe, text = "Select Sites from List", 
            command = function() {
                choose.poll()
                test.sites <<- tk_select.list(sort(unique(poll$siteid)), 
                  multiple = TRUE, title = "2008 Sites")
                tcl("focus", "-force", dialog)
            })
        f.but <<- newguiFilename(dframe, text = "Select Sites from CSV", 
            filter = "{{Comma Separated Values file} {.csv}}")
        tkgrid(s.but, tklabel(dframe, text = "  OR  "), f.but$guiObject)
        tkgrid.configure(s.but, sticky = "nes")
        tkgrid.configure(f.but$guiObject, sticky = "nes")
        tkgrid(tklabel(dframe, text = " "))
        dframe <- guiFrame(borderwidth = 0, sframe = dialog, 
            , relief = "flat")
        tkgrid.configure(dframe, sticky = "news")
        d.but <- tkbutton(dframe, text = "Select directory to save files", 
            command = function() {
                dname <<- tclvalue(tkchooseDirectory())
                tclvalue(lbltxt) <- paste("Current Directory *: ", 
                  dname, sep = "")
            })
        lbltxt <- tclVar("Current Directory *:")
        lbl <- tklabel(dframe, justify = "left", padx = 0, text = tclvalue(lbltxt))
        tkconfigure(lbl, textvariable = lbltxt)
        tkgrid(d.but, lbl)
        tkgrid.configure(d.but, sticky = "nes")
        tkgrid.configure(lbl, sticky = "nes")
        fname <<- guiTextEntry(sframe = dialog, text = "File name *", 
            default = "")
        tkgrid(fname$guiObject)
        dframe <- guiFrame(borderwidth = 0, sframe = dialog, 
            , relief = "flat")
        tkgrid.configure(dframe, sticky = "ns")
        ok.but <- tkbutton(dframe, text = "OK", command = function() {
            data.test("main3")
            if (done == 1) {
                tkdestroy(dialog)
                calc.rembias()
            }
        })
        cancel.but <- tkbutton(dframe, text = "Cancel", command = function() {
            done <<- 3
            flag <- 1
            tkdestroy(dialog)
            stop("cancelled by user")
        })
        tkgrid(ok.but, cancel.but)
        tkgrid.configure(ok.but, sticky = "nes")
        tkgrid.configure(cancel.but, sticky = "nes")
        tcl("focus", "-force", dialog)
        tkwait.window(dialog)
        if (done != 3) {
            test.sites <<- ""
            decide <- tktoplevel()
            tkwm.title(decide, "Decision Time!")
            make.more <- tkbutton(decide, text = "Remove More Sites\nPlots", 
                command = function() {
                  done <<- 0
                  tkdestroy(decide)
                })
            all.done <- tkbutton(decide, text = "All Done", command = function() {
                done <<- 3
                tkdestroy(decide)
            })
            tkgrid(tklabel(decide, text = "Your files are saved.\nDo you want to remove more sites?"))
            tkgrid(make.more, all.done)
            tkgrid.configure(all.done, sticky = "nw")
            tkwait.window(decide)
        }
    }
    rm(list = ls())
    gc()
}