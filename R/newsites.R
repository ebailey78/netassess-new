newsites <- function() 
{
    options(warn = -1)

    pb <- winProgressBar(title = "New Sites Tool Starting Up", 
        min = 0, max = 1, width = 500, label = "Loading data. Please wait")
    for (i in j <- c("o3.pair.summ", "pm25.3day.pair.summ", "pm25.6day.pair.summ", 
        "contpm25.pair.summ", "latlongs", "ozone.dvs", "pm25.dvd", 
        "pm25.dva", "contpm25.dvd", "contpm25.dva")) {
        data(list = i, package = "netassess", envir = environment())
        setWinProgressBar(pb, grep(i, j)/length(j))
    }
    Sys.sleep(2)
    close(pb)
    xypos <- function(xpct, ypct) {
        x <- grconvertX(xpct/100, from = "ndc", to = "user")
        y <- grconvertY(ypct/100, from = "ndc", to = "user")
        return(list(x = x, y = y))
    }
    tolambert <- function(xycoords) {
        tempxy <- SpatialPoints(xycoords, proj4string = CRS("+proj=longlat +ellps=clrk66"))
        tempxy <- spTransform(tempxy, CRS("+proj=lcc +lon_0=97w +lat_0=40n +lat_1=33n +lat_2=45n +ellps=clrk66"))
        dimnames(tempxy@coords)[[2]] <- c("lamx", "lamy")
        return(tempxy@coords)
    }
    tolonglat <- function(xycoords) {
        tempxy <- SpatialPoints(xycoords, proj4string = CRS("+proj=lcc +lon_0=97w +lat_0=40n +lat_1=33n +lat_2=45n +ellps=clrk66"))
        tempxy <- spTransform(tempxy, CRS("+proj=longlat +ellps=clrk66"))
        dimnames(tempxy@coords)[[2]] <- c("longitude", "latitude")
        return(tempxy@coords)
    }
    makesites <- function(maininput) {
        if (grepl("ozone", maininput)) {
            dvs <<- ozone.dvs
            pair.summ <- o3.pair.summ
            parm <<- "Ozone"
        }
        if (grepl("pm25frm", maininput)) {
            parm <<- "PM2.5 (FRM)"
            if (grepl("3", maininput)) {
                pair.summ <- pm25.3day.pair.summ
                sampsched <<- "3-day"
            }
            else {
                pair.summ <- pm25.6day.pair.summ
                sampsched <<- "6-day"
            }
            if (grepl("annual", maininput)) {
                dvs <<- pm25.dva
                pmstd <<- "Annual"
            }
            else {
                dvs <<- pm25.dvd
                pmstd <<- "Daily"
            }
        }
        if (grepl("contpm25", maininput)) {
            parm <<- "Continuous PM2.5"
            pair.summ <- contpm25.pair.summ
            if (grepl("annual", maininput)) {
                dvs <<- contpm25.dva
                pmstd <<- "Annual"
            }
            else {
                dvs <<- contpm25.dvd
                pmstd <<- "Daily"
            }
        }
        pair.summ$sitecomb <- sapply(strsplit(pair.summ$sitecomb, 
            split = " ", fixed = TRUE), function(x) paste(sort(x), 
            collapse = " "))
        guicorr <- as.numeric(tclvalue(pair.corr$object))
        guidist <- as.numeric(tclvalue(distance$object))
        guigrad <- as.numeric(tclvalue(gradient))
        sig.pair.summ <- pair.summ[pair.summ$pvalue <= 0.025 & 
            pair.summ$corr^2 <= guicorr & pair.summ$dist >= guidist & 
            pair.summ$dist <= 500 * 1.609344, ]
        sig.pair.summ <- sig.pair.summ[is.na(sig.pair.summ$site1) == 
            FALSE, ]
        sig.pair.summ <- sig.pair.summ[abs(sig.pair.summ$Mean) >= 
            guigrad, ]
        if (nrow(sig.pair.summ) == 0) {
            tkmessageBox(title = "Gradient Out of Bounds!", message = paste("The gradient value is out of bounds!\nPlease try a value between 0 and ", 
                uprbnd, sep = ""), icon = "info", type = "ok")
            stop
        }
        if (parm == "Ozone") {
            uprbnd <- round(max(abs(sig.pair.summ$Mean), na.rm = TRUE), 
                digits = 3)
        }
        else {
            uprbnd <- round(max(abs(sig.pair.summ$Mean), na.rm = TRUE), 
                digits = 1)
        }
        sig.pair.summ$place.dist <- sqrt((sig.pair.summ$lamx1 - 
            sig.pair.summ$lamx2)^2 + (sig.pair.summ$lamy1 - sig.pair.summ$lamy2)^2)/2
        m <- (sig.pair.summ$lamy2 - sig.pair.summ$lamy1)/(sig.pair.summ$lamx2 - 
            sig.pair.summ$lamx1)
        int <- sig.pair.summ$lamy2 - (sig.pair.summ$lamx2 * m)
        a <- 1 + m^2
        b <- (2 * m * int) - (2 * sig.pair.summ$lamx1) - (2 * 
            sig.pair.summ$lamy1 * m)
        c <- sig.pair.summ$lamx1^2 + sig.pair.summ$lamy1^2 - 
            (2 * sig.pair.summ$lamy1 * int) + int^2 - sig.pair.summ$place.dist^2
        angle <- atan2(sig.pair.summ$lamy2 - sig.pair.summ$lamy1, 
            sig.pair.summ$lamx2 - sig.pair.summ$lamx1)
        angle <- ifelse(angle < 0, angle + (2 * pi), angle)
        sig.pair.summ$new.x <- ifelse(pi/2 <= angle & angle <= 
            3 * pi/2, (-b - sqrt(b^2 - (4 * a * c)))/(2 * a), 
            (-b + sqrt(b^2 - (4 * a * c)))/(2 * a))
        sig.pair.summ$new.y <- int + (m * sig.pair.summ$new.x)
        sig.pair.summ <- cbind(sig.pair.summ, tolonglat(sig.pair.summ[, 
            c("new.x", "new.y")]))
        sig.new.sites <- unique(sig.pair.summ[, c("sitecomb", 
            "longitude", "latitude")])
        sig.new.sites <- cbind(sig.new.sites, tolambert(unique(sig.pair.summ[, 
            c("longitude", "latitude")])))
        the.grid <- GridTopology(c(-2406000, -1614000), rep(12000, 
            2), c(396, 246))
        the.grid <- SpatialGrid(the.grid, proj4string = CRS("+proj=lcc +lon_0=97w +lat_0=40n +lat_1=33n +lat_2=45n"))
        temp <- SpatialPoints(sig.new.sites[, c("lamx", "lamy")], 
            proj4string = CRS("+proj=lcc +lon_0=97w +lat_0=40n +lat_1=33n +lat_2=45n"))
        sig.new.sites <- SpatialPointsDataFrame(temp, sig.new.sites)
        sig.new.sites <- data.frame(coordinates(the.grid)[overlay(the.grid, 
            sig.new.sites), ], sig.new.sites@data)
        sig.new.sites <- merge(dvs, sig.new.sites, by.x = c("cell.x", 
            "cell.y"), by.y = c("s1", "s2"))
        sig.new.sites <- sig.new.sites[ifelse(is.na(sig.new.sites$instates), 
            NA, sig.new.sites$prob) >= as.numeric(tclvalue(probcut$object))/100, 
            c("sitecomb", "lamx", "lamy")]
        temp <- sig.new.sites$sitecomb
        sig.new.sites <<- sig.new.sites[, c("lamx", "lamy")]
        non.sig.sites <- pair.summ[!pair.summ$sitecomb %in% sig.new.sites$sitecomb, 
            c("site1", "site2")]
        non.sig.sites <<- list(non.sig.sites = unique(c(non.sig.sites$site1, 
            non.sig.sites$site2)), sig.sites = unique(unlist(strsplit(temp, 
            split = " ", fixed = TRUE))))
        site1 <- unique(pair.summ[, c("site1", "long1", "lat1")])
        names(site1) <- c("site", "long", "lat")
        site2 <- unique(pair.summ[, c("site2", "long2", "lat2")])
        names(site2) <- c("site", "long", "lat")
        ozone.sites <<- unique(rbind(site1, site2))
        rm(list = c("site1", "site2"))
        gc()
        try(tkrreplot(img), silent = TRUE)
    }
    makemap <- function() {
        ozone.dvs0608 <- matrix(ifelse(is.na(dvs$instates), NA, 
            dvs$prob), ncol = max(dvs$row))
        breakpts <- seq(0, 1, length.out = 1001)
        col.palette <- colorRampPalette(c(rgb(1, 1, 1), rgb(1, 
            1, 0), "#FF8C00", rgb(1, 0, 0), rgb(0.8, 0, 0.8)))(1000)
        par(mar = c(5, 0, 0, 0), bg = "white", xpd = TRUE)
        states <- map2SpatialLines(map("state", plot = FALSE), 
            proj4string = CRS("+proj=longlat +ellps=clrk66"))
        states <- spTransform(states, CRS("+proj=lcc +lon_0=97w +lat_0=40n +lat_1=33n +lat_2=45n"))
        image(unique(dvs$cell.x), unique(dvs$cell.y), ozone.dvs0608, 
            xlab = "", ylab = "", axes = FALSE, breaks = breakpts, 
            col = col.palette)
        plot(states, col = "gray", add = TRUE)
        gradient.rect(xypos(5, 5)$x, xypos(5, 5)$y, xypos(95, 
            10)$x, xypos(95, 10)$y, col = col.palette)
        text(xypos(seq(5, 95, length.out = 11), rep(4, 11)), 
            labels = seq(0, 100, length.out = 11))
        old.sites1 <- ozone.sites[ozone.sites$site %in% non.sig.sites$non.sig.sites, 
            ]
        old.sites1 <- tolambert(old.sites1[, c("long", "lat")])
        old.sites2 <- ozone.sites[ozone.sites$site %in% non.sig.sites$sig.sites, 
            ]
        old.sites2 <- tolambert(old.sites2[, c("long", "lat")])
        points(old.sites1, pch = 19, col = "gray")
        points(old.sites2, pch = 19, col = "black")
        points(sig.new.sites, pch = 2, col = "blue")
        title.text <- paste("2008", parm, "sites", sep = " ")
        if (parm == "PM2.5 (FRM)") {
            title.text <- paste(title.text, " (", sampsched, 
                ")", sep = "")
        }
        if (parm %in% c("PM2.5 (FRM)", "Continuous PM2.5")) {
            units <- "ug/m3"
        }
        else {
            units <- "ppm"
        }
        text(x = xypos(50, 97)$x, y = xypos(50, 97)$y, labels = title.text)
        text(x = xypos(60, 95)$x, y = xypos(60, 95)$y, labels = paste("Site Pair Correlation < ", 
            tclvalue(pair.corr$object), sep = ""), pos = 4, cex = 0.85)
        text(x = xypos(60, 93)$x, y = xypos(60, 93)$y, labels = paste("Minimum Distance between Site Pairs: ", 
            tclvalue(distance$object), "km", sep = ""), pos = 4, 
            cex = 0.85)
        text(x = xypos(60, 91)$x, y = xypos(60, 91)$y, labels = paste("Difference between Site Pairs: ", 
            tclvalue(gradient), " ", units, sep = ""), pos = 4, 
            cex = 0.85)
        text(x = xypos(60, 89)$x, y = xypos(60, 89)$y, labels = paste("Probability of Exceeding 85% of NAAQS: ", 
            tclvalue(probcut$object), "%", sep = ""), pos = 4, 
            cex = 0.85)
        points(x = xypos(5, 21)$x, y = xypos(5, 21)$y, pch = 19, 
            col = "black")
        text(x = xypos(5, 21)$x, y = xypos(5, 21)$y, labels = "Existing Sites from Site Pairs Meeting Criteria", 
            pos = 4)
        points(x = xypos(5, 18)$x, y = xypos(5, 18)$y, pch = 19, 
            col = "gray")
        text(x = xypos(5, 18)$x, y = xypos(5, 18)$y, labels = "Existing Sites from Site Pairs Not Meeting Criteria", 
            pos = 4)
        points(x = xypos(5, 15)$x, y = xypos(5, 15)$y, pch = 2, 
            col = "blue")
        text(x = xypos(5, 15)$x, y = xypos(5, 15)$y, labels = "Possible New Sites", 
            pos = 4)
        text(x = xypos(50, 12)$x, y = xypos(50, 12)$y, labels = "Probability of Exceeding 85% of NAAQS")
    }
    data.test <- function(checkvar) {
        done <<- 0
        if (checkvar == "main3") {
            if (tclvalue(fname$object) == "" | dname == "" | 
                tclvalue(tcl(polltree, "selection", "get")) == 
                  "") {
                tt <- tkmessageBox(title = "'Doh!'", message = "One or more of the items marked by '*' is missing!", 
                  icon = "info", type = "ok")
            }
            else {
                done <<- 1
            }
        }
    }
    savefiles <- function() {
        pdf(file = paste(dname, "/", tclvalue(fname$object), 
            ".pdf", sep = ""), width = 12, height = 10.5)
        makemap()
        dev.off()
        kml <- file(paste(dname, "/", tclvalue(fname$object), 
            ".kml", sep = ""), open = "w")
        cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", "<kml xmlns=\"http://earth.google.com/kml/2.1\">", 
            "<Document>", paste("\t<name>Possible ", parm, " \"New\" Sites</name>", 
                sep = ""), "\t<visibility>0</visibility>", "\t<open>0</open>", 
            "\t<Style id=\"yessites\">", "\t\t<IconStyle>", "\t\t\t<scale>1</scale>", 
            "\t\t\t<color>ff0000ff</color>", "\t\t\t<Icon>", 
            "\t\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>", 
            "\t\t\t</Icon>", "\t\t</IconStyle>", "\t</Style>", 
            "\t<Style id=\"nosites\">", "\t\t<IconStyle>", "\t\t\t<scale>1</scale>", 
            "\t\t\t<color>ff808080</color>", "\t\t\t<Icon>", 
            "\t\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>", 
            "\t\t\t</Icon>", "\t\t</IconStyle>", "\t</Style>", 
            "\t<Style id=\"newsites\">", "\t\t<IconStyle>", "\t\t\t<scale>1</scale>", 
            "\t\t\t<color>ffff0000</color>", "\t\t\t<Icon>", 
            "\t\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/star.png</href>", 
            "\t\t\t</Icon>", "\t\t</IconStyle>", "\t</Style>", 
            "\t<Style>", "\t\t<LineStyle id=\"connect\">", "\t\t\t<color>ffffffff</color>", 
            "\t\t\t<colorMode>normal</colorMode>", "\t\t\t<width>1</width>", 
            "\t\t</LineStyle>", "\t</Style>", "\t<Style id=\"noitems\">", 
            "\t\t<ListStyle>", "\t\t\t<listItemType>checkHideChildren</listItemType>", 
            "\t\t</ListStyle>", "\t</Style>", "\t<Folder>", "\t\t<name>Existing Ozone Sites</name>", 
            "\t\t<visibility>0</visibility>", "\t\t<open>0</open>", 
            "\t\t<styleUrl>#noitems</styleUrl>\n", sep = "\n", 
            file = kml)
        yessites <- ozone.sites[ozone.sites$site %in% non.sig.sites$sig.sites, 
            ]
        nosites <- ozone.sites[ozone.sites$site %in% non.sig.sites$non.sig.sites[!non.sig.sites$non.sig.sites %in% 
            non.sig.sites$sig.sites], ]
        cat(paste("\t\t<Placemark>\n\t\t\t<styleUrl>#yessites</styleUrl>\n\t\t\t<Point>\n\t\t\t\t<extrude>0</extrude>\n\t\t\t\t<altitudeMode>ClampToGround</altitudeMode>\n\t\t\t\t<coordinates>", 
            yessites$long, ",", yessites$lat, ",0</coordinates>\n\t\t\t</Point>\n\t\t</Placemark>\n", 
            sep = ""), paste("\t\t<Placemark>\n\t\t\t<styleUrl>#nosites</styleUrl>\n\t\t\t<Point>\n\t\t\t\t<extrude>0</extrude>\n\t\t\t\t<altitudeMode>ClampToGround</altitudeMode>\n\t\t\t\t<coordinates>", 
            nosites$long, ",", nosites$lat, ",0</coordinates>\n\t\t\t</Point>\n\t\t</Placemark>\n", 
            sep = ""), file = kml)
        cat("\t</Folder>", "\t<Folder>", "\t\t<name>Possible New Sites</name>", 
            "\t\t<visibility>0</visibility>", "\t\t<open>0</open>", 
            "\t\t<styleUrl>#noitems</styleUrl>\n", sep = "\n", 
            file = kml)
        cat(paste("\t\t<Placemark>\n\t\t\t<styleUrl>#newsites</styleUrl>\n\t\t\t<Point>\n\t\t\t\t<extrude>0</extrude>\n\t\t\t\t<altitudeMode>ClampToGround</altitudeMode>\n\t\t\t\t<coordinates>", 
            sig.new.sites$long, ",", sig.new.sites$lat, ",0</coordinates>\n\t\t\t</Point>\n\t\t</Placemark>\n", 
            sep = ""), file = kml)
        cat("\t</Folder>", "</Document>", "</kml>", sep = "\n", 
            file = kml)
        close(con = kml)
        yessites$type <- "meets criteria"
        nosites$type <- "does not meet critera"
        sig.new.sites$type <- "possible new site"
        sig.new.sites$site <- paste("site", 1:nrow(sig.new.sites), 
            sep = "")
        names(sig.new.sites)[names(sig.new.sites) %in% c("longitude", 
            "latitude")] <- c("long", "lat")
        finalsites <- rbind(yessites, nosites, sig.new.sites)
        write.csv(finalsites, file = paste(dname, "/", tclvalue(fname$object), 
            ".csv", sep = ""), row.names = FALSE)
        temp <- SpatialPoints(finalsites[, c("long", "lat")], 
            proj4string = CRS("+proj=longlat +ellps=clrk66"))
        finalsites <- SpatialPointsDataFrame(temp, finalsites)
        writePointsShape(finalsites, fn = paste(dname, "/", tclvalue(fname$object), 
            sep = ""))
        tkmessageBox(title = "", message = "All the files are saved.", 
            icon = "info", type = "ok")
    }
    varin <<- "ozone8hr"
    dname <- ""
    dialog <- tktoplevel()
    tclRequire("BWidget")
    sw <- tclvalue(tkwinfo("screenwidth", dialog))
    sh <- tclvalue(tkwinfo("screenheight", dialog))
    tkwm.maxsize(dialog, sw, sh)
    tkwm.state(dialog, "zoomed")
    tkwm.title(dialog, "Possible Locations for New Sites")
    swin <- tkwidget(dialog, "ScrolledWindow")
    tkpack(swin, fill = "both", expand = "yes")
    c1 <- tkcanvas(swin, scrollregion = c(0, 0, 2000, 2000))
    tcl(swin, "setwidget", c1)
    fr <- tkframe(c1)
    tkcreate(c1, "window", 0, 0, anchor = "nw", window = fr)
    fr1 <- tkframe(fr)
    fr2 <- tkframe(fr)
    tkgrid(fr1, fr2, sticky = "news")
    reset.but <- tkbutton(fr1, text = "Reset", command = function() {
        tclvalue(pair.corr$object) <<- 0.5
        tclvalue(distance$object) <<- 100
        tclvalue(gradient) <<- 0
        tclvalue(probcut$object) <<- 80
        makesites(varin)
    })
    guiSet("SLIDER_LENGTH", 100)
    tkgrid(tklabel(fr1, text = "Valid choices in the tree below are colored \"blue\".\nYou may have to expand an individual branch\nof the tree by clicking on the \"plus\" to\nget to a valid choice. *", 
        justify = "left"), sticky = "nw")
    dframe <- guiFrame(fr1)
    tkgrid(dframe, sticky = "nw")
    scrx <- tkscrollbar(dframe, command = function(...) tkxview(polltree, 
        ...), orient = "horizontal")
    scry <- tkscrollbar(dframe, command = function(...) tkyview(polltree, 
        ...), orient = "vertical")
    polltree <- tkwidget(dframe, "Tree", xscrollcommand = function(...) tkset(scrx, 
        ...), yscrollcommand = function(...) tkset(scry, ...))
    tkgrid(polltree, scry)
    tkgrid.configure(scry, sticky = "news")
    tkgrid(scrx, sticky = "ew")
    tkinsert(polltree, "end", "root", "ozone", text = "Ozone", 
        selectable = "no")
    tkinsert(polltree, "end", "ozone", "ozone8hr", text = "8 Hour Std", 
        fill = "blue")
    tkinsert(polltree, "end", "root", "pm25frm", text = "PM2.5 (FRM)", 
        selectable = "no")
    tkinsert(polltree, "end", "pm25frm", "day3", text = "3-day schedule", 
        selectable = "no")
    tkinsert(polltree, "end", "pm25frm", "day6", text = "6-day schedule", 
        selectable = "no")
    tkinsert(polltree, "end", "day3", "pm25frm3annual", text = "Annual Std", 
        fill = "blue")
    tkinsert(polltree, "end", "day3", "pm25frm3daily", text = "Daily Std", 
        fill = "blue")
    tkinsert(polltree, "end", "day6", "pm25frm6annual", text = "Annual Std", 
        fill = "blue")
    tkinsert(polltree, "end", "day6", "pm25frm6daily", text = "Daily Std", 
        fill = "blue")
    tkinsert(polltree, "end", "root", "contpm25", text = "Continuous PM2.5", 
        selectable = "no")
    tkinsert(polltree, "end", "contpm25", "contpm25annual", text = "Annual Std", 
        fill = "blue")
    tkinsert(polltree, "end", "contpm25", "contpm25daily", text = "Daily Std", 
        fill = "blue")
    tcl(polltree, "bindText", "<ButtonRelease-1>", function(...) if (tclvalue(tcl(polltree, 
        "selection", "get")) != "") {
        varin <<- tclvalue(tcl(polltree, "selection", "get"))
        makesites(varin)
    })
    pair.corr <<- guiSlider(sframe = fr1, text = "Correlation between Site Pairs (R2)", 
        default = 0.5, min = 0, max = 1, step = 0.01, update = function(...) makesites(varin, 
            ...))
    tkgrid(pair.corr$guiObject, reset.but, sticky = "nws")
    distance <<- guiSlider(sframe = fr1, text = "Minimum Distance between Site Pairs (km)", 
        default = 100, min = 0, max = 500, step = 50, update = function(...) makesites(varin, 
            ...))
    tkgrid(distance$guiObject, sticky = "nw")
    dframe <- guiFrame(sframe = fr1)
    tkgrid.configure(dframe, sticky = "nesw")
    gradient <<- tclVar("0")
    gradlbl <- tkentry(dframe, width = 10, textvariable = gradient)
    gradlbltxt <- tklabel(dframe, text = "Average Difference(Gradient) between Site Pairs")
    tkgrid(gradlbltxt, gradlbl)
    tkgrid.configure(gradlbltxt, sticky = "nes")
    tkgrid.configure(gradlbl, sticky = "nes")
    tkbind(gradlbl, "<Return>", function(...) makesites(varin, 
        ...))
    probcut <<- guiSlider(sframe = fr1, text = "Probability of Design Value Exceeding 85% of the NAAQS", 
        default = 80, min = 0, max = 100, step = 10, update = function(...) makesites(varin, 
            ...))
    tkgrid(probcut$guiObject)
    tkgrid.configure(probcut$guiObject, sticky = "nw")
    dframe <- guiFrame(sframe = fr2)
    tkgrid(dframe, sticky = "news")
    d.but <- tkbutton(dframe, text = "Select directory to save files", 
        command = function() {
            dname <<- tclvalue(tkchooseDirectory())
            tclvalue(lbltxt) <- paste("Current Directory *: ", 
                dname, sep = "")
        })
    lbltxt <- tclVar("Current Directory *:")
    lbl <- tklabel(dframe, text = tclvalue(lbltxt))
    tkconfigure(lbl, textvariable = lbltxt)
    tkgrid(d.but, lbl)
    tkgrid.configure(d.but, sticky = "nes")
    tkgrid.configure(lbl, sticky = "nes")
    fname <<- guiTextEntry(sframe = fr2, text = "File name *", 
        default = "")
    tkgrid(fname$guiObject)
    tkgrid.configure(fname$guiObject, sticky = "nesw")
    makesites(varin)
    img <<- tkrplot(fr2, fun = makemap, hscale = 2.3333, vscale = 2)
    tkgrid(img)
    dframe2 <- guiFrame(fr2)
    tkgrid(dframe2, sticky = "news")
    save.but <- tkbutton(dframe2, text = "Save plot", command = function() {
        data.test("main3")
        if (done == 1) {
            savefiles()
        }
    })
    cancel.but <- tkbutton(dframe2, text = "Cancel", command = function() {
        done <<- 2
        tkdestroy(dialog)
        stop("cancelled by user")
    })
    done.but <- tkbutton(dframe2, text = "Done", command = function() {
        done <<- 2
        tkdestroy(dialog)
    })
    tkgrid(save.but, done.but, cancel.but)
    tkgrid.configure(cancel.but, sticky = "nw")
    tkconfigure(c1, scrollregion = c(0, 0, as.numeric(tkwinfo("width", 
        fr)), as.numeric(tkwinfo("height", fr))))
    tcl("focus", "-force", dialog)
    tkwait.window(dialog)
    rm(list = ls(pos = 1), pos = 1)
    gc(verbose = FALSE)
}