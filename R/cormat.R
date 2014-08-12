cormat <- function() 
{
    options(warn = -1)
    require(maps)
    require(sp)
    require(rgdal)
    require(maptools)
    require(tcltk)
    require(tkrplot)
    require(fgui)
    require(sqldf)
    require(Hmisc)
    require(ellipse)
    require(plotrix)
    pb <- winProgressBar(title = "Correlation Matrix Tool Starting Up", 
        min = 0, max = 1, width = 500, label = "Loading data. Please wait")
    for (i in j <- c("ozone", "pm25avg", "contpm25_all", "latlongs")) {
        data(list = i, package = "netassess")
        setWinProgressBar(pb, grep(i, j)/length(j))
    }
    Sys.sleep(2)
    close(pb)
    counties <- map("county", plot = FALSE)
    states <- map("state", plot = FALSE)
    tolambert <- function(xycoords) {
        tempxy <- SpatialPoints(xycoords, proj4string = CRS("+proj=longlat +ellps=clrk66"))
        tempxy <- spTransform(tempxy, CRS("+proj=lcc +lon_0=97w +lat_0=40n +lat_1=33n +lat_2=45n +ellps=clrk66"))
        dimnames(tempxy@coords)[[2]] <- c("lamx", "lamy")
        return(tempxy@coords)
    }
    xypos <- function(xpct, ypct) {
        x <- grconvertX(xpct/100, from = "ndc", to = "user")
        y <- grconvertY(ypct/100, from = "ndc", to = "user")
        return(list(x = x, y = y))
    }
    makestates <- function() {
        par(mar = c(0, 0, 3, 0), bg = "white")
        if (length(boxcoords) == 2) {
            temp <- matrix(unlist(imgbox), ncol = 2)
            mapbox <<- t(apply(temp, 1, sort))
        }
        map("county", col = "gray", xlim = mapbox[1, ], ylim = mapbox[2, 
            ])
        lines(states, col = "black")
        if (length(boxcoords) == 0) {
            mtext("Click on the map to position the lower left corner of the box", 
                col = "blue")
        }
        if (length(boxcoords) == 1) {
            mtext("Click on the map to position the upper right corner of the box", 
                col = "red")
            points(imgbox[[1]], col = "blue", pch = 15, cex = 1.5)
        }
        if (length(boxcoords) == 2) {
            rect(imgbox[[1]][1], imgbox[[1]][2], imgbox[[2]][1], 
                imgbox[[2]][2], border = "red")
            mapbox
        }
        plotsize <<- par("plt")
        usrcoords <<- matrix(par("usr"), ncol = 2, byrow = TRUE)
    }
    makeclick <- function(x, y) {
        if (length(boxcoords) == 2) {
            boxcoords <<- list()
            imgbox <<- list()
        }
        height <- as.numeric(tkwinfo("reqheight", img))
        width <- as.numeric(tkwinfo("reqwidth", img))
        boxcoords[[length(boxcoords) + 1]] <<- c(x, height - 
            as.numeric(y))
        extent <- matrix(plotsize * rep(c(width, height), each = 2), 
            ncol = 2, byrow = TRUE)
        rtotc <- apply(usrcoords, 1, function(x) x[2] - x[1])/apply(extent, 
            1, function(x) x[2] - x[1])
        imgbox <<- lapply(boxcoords, function(z) {
            matrix((as.numeric(z) - extent[, 1]) * rtotc + usrcoords[, 
                1], ncol = 2)
        })
        if (length(boxcoords) >= 1) {
            tkrreplot(img)
        }
    }
    choosearea <- function(mons) {
        boxcoords <<- list()
        imgbox <<- list()
        mapbox <<- matrix(c(min(counties$x, na.rm = TRUE), max(counties$x, 
            na.rm = TRUE), min(counties$y, na.rm = TRUE), max(counties$y, 
            na.rm = TRUE)), ncol = 2, byrow = TRUE)
        tt <- tktoplevel()
        tcl("focus", "-force", tt)
        tkwm.title(tt, "Choose an Area")
        img <<- tkrplot(tt, fun = makestates, hscale = 2, vscale = 2)
        tkpack(img)
        resetmap <- function() {
            boxcoords <<- list()
            imgbox <<- list()
            mapbox <<- matrix(c(min(counties$x, na.rm = TRUE), 
                max(counties$x, na.rm = TRUE), min(counties$y, 
                  na.rm = TRUE), max(counties$y, na.rm = TRUE)), 
                ncol = 2, byrow = TRUE)
            tkrreplot(img)
        }
        tkbind(img, "<Button-1>", makeclick)
        topmenu <- tkmenu(tt)
        tkconfigure(tt, menu = topmenu)
        tkadd(topmenu, "command", label = "Reset Map", command = resetmap)
        tkadd(topmenu, "command", label = "Done Zooming", command = function() {
            data.test("imgbox")
            if (done == 1) {
                tkdestroy(tt)
            }
        })
        tcl("focus", "-force", tt)
        tkwait.window(tt)
    }
    data.test <- function(checkvar) {
        done <<- 0
        if (checkvar == "imgbox") {
            if (length(imgbox) == 0) {
                tt <- tkmessageBox(title = "'Doh!'", message = "You must choose an area from the map!", 
                  icon = "info", type = "ok")
            }
            else {
                done <<- 1
            }
        }
        if (checkvar == "main3") {
            if (tclvalue(yr$object) == "" | tclvalue(fname$object) == 
                "" | dname == "") {
                tt <- tkmessageBox(title = "'Doh!'", message = "One or more of the items marked by '*' is missing!", 
                  icon = "info", type = "ok")
                tcl("focus", "-force", dialog)
            }
            else {
                done <<- 1
            }
        }
    }
    make.cormat <- function() {
        pb <- winProgressBar(title = "Calculating", label = "0% done", 
            max = 100)
        temp <- matrix(unlist(imgbox), ncol = 2)
        temp <- t(apply(temp, 1, sort))
        area_sites <- latlongs[temp[1, 1] <= latlongs$longitude & 
            latlongs$longitude <= temp[1, 2] & temp[2, 1] <= 
            latlongs$latitude & latlongs$latitude <= temp[2, 
            2], ]
        if (grepl("ozone", varin)) {
            crit.val <- "maxval"
            poll <- ozone[as.numeric(substr(ozone$colldate, 6, 
                7)) %in% 5:9 & ozone$yr %in% eval(parse(text = tclvalue(yr$object))) & 
                ozone$siteid %in% area_sites$siteid, ]
            if (length(unique(poll$siteid)) == 0) {
                close(pb)
                return(0)
            }
            compsites <- sqldf(paste("select distinct siteid\n\t\t\t\t\t\t\t\t\tfrom (select distinct siteid, count(*) as nyears \n\t\t\t\t\t\t\t\t\t\t\tfrom (select distinct yr, siteid, count(*) as ndays \n\t\t\t\t\t\t\t\t\t\t\t\t\tfrom poll where (obspct>=75 or maxval>0.075) and yr in (", 
                paste(eval(parse(text = tclvalue(yr$object))), 
                  collapse = ","), ") group by yr, siteid) \n\t\t\t\t\t\t\t\t\t\t\twhere ndays>=153*0.75 \n\t\t\t\t\t\t\t\t\t\t\tgroup by siteid)\n\t\t\t\t\t\t\t\t\twhere nyears=", 
                length(eval(parse(text = tclvalue(yr$object)))), 
                sep = ""), method = "raw")$siteid
        }
        else if (grepl("pm25frm|contpm25", varin)) {
            poll <- pm25avg
            crit.val <- "pm25"
            startday <- list("2005-01-01", "2005-01-04")
            names(startday) <- c("3-day", "6-day")
            if (grepl("3-day", varin)) {
                sampsched <- "3-day"
            }
            else {
                sampsched <- "6-day"
            }
            daysuse <- as.character(seq(from = as.Date(startday[[sampsched]], 
                "%Y-%m-%d"), to = as.Date("2008-12-31", "%Y-%m-%d"), 
                by = as.numeric(substr(sampsched, 1, 1))))
            if (grepl("contpm25", varin)) {
                daysuse <- as.character(seq(from = as.Date("2005-01-01", 
                  "%Y-%m-%d"), to = as.Date("2008-12-31", "%Y-%m-%d"), 
                  by = 1))
                poll <- contpm25_all
                crit.val <- "contpm25"
            }
            poll <- poll[poll$colldate %in% daysuse & poll$yr %in% 
                eval(parse(text = tclvalue(yr$object))) & poll$siteid %in% 
                area_sites$siteid, ]
            poll$qtr <- ceiling(as.numeric(substr(poll$colldate, 
                6, 7))/3)
            if (length(unique(poll$siteid)) == 0) {
                close(pb)
                return(0)
            }
            daysuse <- daysuse[as.numeric(substr(daysuse, 1, 
                4)) %in% eval(parse(text = tclvalue(yr$object)))]
            daysuse <- data.frame(days = daysuse, yr = as.numeric(substr(daysuse, 
                1, 4)), qtr = ceiling(as.numeric(substr(daysuse, 
                6, 7))/3), stringsAsFactors = FALSE)
            ndays <- sqldf("select distinct yr, qtr, count(*) as ndays\n\t\t\t\t\tfrom daysuse\n\t\t\t\t\tgroup by yr, qtr", 
                method = "raw")
            compsites <- sqldf("select distinct siteid, yr, qtr, count(*) as nobs\n\t\t\t\t\tfrom poll\n\t\t\t\t\tgroup by siteid, yr, qtr", 
                method = "raw")
            compsites <- merge(compsites, ndays, by = c("yr", 
                "qtr"), all.x = TRUE)
            compsites <- sqldf(paste("select distinct siteid from \n\t\t\t\t\t\t\t(select distinct siteid, count(*) as nyrs \n\t\t\t\t\t\t\t\tfrom (select distinct siteid, yr, count(*) as nqtrs\n\t\t\t\t\t\t\t\t\t\tfrom compsites\n\t\t\t\t\t\t\t\t\t\twhere nobs>=0.75*ndays\n\t\t\t\t\t\t\t\t\t\tgroup by siteid, yr)\n\t\t\t\t\t\t\t\twhere nqtrs=4\n\t\t\t\t\t\t\t\tgroup by siteid)\n\t\t\t\t\t\twhere nyrs=", 
                length(eval(parse(text = tclvalue(yr$object)))), 
                sep = ""), method = "raw")$siteid
        }
        setWinProgressBar(pb, value = 10, title = "Calculating", 
            label = "10% done")
        if (length(compsites) == 0) {
            close(pb)
            return(0)
        }
        setWinProgressBar(pb, value = 20, title = "Calculating", 
            label = "20% done")
        area_poll <- poll[poll$siteid %in% compsites, c("siteid", 
            "colldate", crit.val)]
        test <- reshape(area_poll, idvar = "colldate", timevar = "siteid", 
            direction = "wide")
        test2 <- try(cor(test[, grep(crit.val, names(test))], 
            method = "pearson", use = "pairwise.complete.obs"), 
            silent = TRUE)
        setWinProgressBar(pb, value = 50, title = "Calculating", 
            label = "50% done")
        if (!inherits(test2, "try-error")) {
            colnames(test2) <- do.call("rbind", strsplit(colnames(test2), 
                ".", fixed = TRUE))[, 2]
            rownames(test2) <- do.call("rbind", strsplit(rownames(test2), 
                ".", fixed = TRUE))[, 2]
            varpairs <- as.matrix(expand.grid(site1 = unique(area_poll$siteid), 
                site2 = unique(area_poll$siteid)))
            area_poll_sites <- latlongs[latlongs$siteid %in% 
                unique(area_poll$siteid), ]
            lamb_coord <- tolambert(area_poll_sites[, c("longitude", 
                "latitude")])
            area_poll_sites <- cbind(area_poll_sites, lamb_coord)
            rat_stats <- apply(varpairs, 1, function(z) {
                val1 <- eval(parse(text = paste(paste("test$", 
                  crit.val, sep = ""), z[1], sep = ".")))
                val2 <- eval(parse(text = paste(paste("test$", 
                  crit.val, sep = ""), z[2], sep = ".")))
                rats <- abs(val2 - val1)/mean(c(val2, val1), 
                  na.rm = TRUE)
                ratio <- mean(rats, na.rm = TRUE)
                med_rat <- median(rats, na.rm = TRUE)
                sd_rat <- sd(rats, na.rm = TRUE)
                min_rat <- min(rats, na.rm = TRUE)
                max_rat <- max(rats, na.rm = TRUE)
                nobs <- length(rats[!is.na(rats)])
                return(data.frame(site1 = z[1], site2 = z[2], 
                  ratio, med_rat, sd_rat, min_rat, max_rat, nobs, 
                  stringsAsFactors = FALSE))
            })
            rat_stats <- do.call("rbind", rat_stats)
            o3_ratio <- rat_stats$ratio
            colors <- colorRampPalette(c("white", "yellow", "orange", 
                "red", "purple"))(1000)
            colors <- c(colors, colors[1000])
            breakpts <- cut(o3_ratio, c(seq(0, 1, length.out = 1001), 
                9999), labels = colors, include.lowest = TRUE, 
                right = TRUE)
            breakpts <- matrix(breakpts, ncol = nrow(area_poll_sites), 
                nrow = nrow(area_poll_sites), dimnames = rep(list(unique(area_poll_sites$siteid)), 
                  2))
            o3_distance <- apply(varpairs, 1, function(z) {
                site1 <- area_poll_sites[area_poll_sites$siteid == 
                  z[1], ]
                site2 <- area_poll_sites[area_poll_sites$siteid == 
                  z[2], ]
                distance <- round(sqrt((site1$lamx - site2$lamx)^2 + 
                  (site1$lamy - site2$lamy)^2)/1000, digits = 0)
                return(distance)
            })
            names(rat_stats) <- c("site1", "site2", "avg_rel_diff", 
                "median_rel_diff", "sd_rel_diff", "min_rel_diff", 
                "max_rel_diff", "nobs")
            csv.out <- data.frame(rat_stats, corr = as.vector(test2), 
                distance = o3_distance, stringsAsFactors = FALSE)
            write.csv(csv.out, file = paste(dname, "/", tclvalue(fname$object), 
                ".csv", sep = ""), row.names = FALSE)
            o3_distance <- matrix(as.character(o3_distance), 
                ncol = nrow(area_poll_sites), nrow = nrow(area_poll_sites), 
                dimnames = rep(list(unique(area_poll_sites$siteid)), 
                  2))
            o3_distance[lower.tri(o3_distance)] <- ""
            data.out <- data.frame(site1 = varpairs[, 1], site2 = varpairs[, 
                2], corr = as.vector(test2), avgratio = as.vector(o3_ratio), 
                distance = as.numeric(as.vector(t(o3_distance))), 
                stringsAsFactors = FALSE)
            data.out <- data.out[!is.na(data.out$distance), ]
        }
        else {
            close(pb)
            return(0)
        }
        setWinProgressBar(pb, value = 90, title = "Calculating", 
            label = "90% done")
        pdf(file = paste(dname, "/", tclvalue(fname$object), 
            ".pdf", sep = ""), width = 12, height = 10.5)
        plotcorr(test2, col = breakpts, type = "lower", diag = TRUE, 
            cex.lab = 0.75, mar = c(1, 5, 0, 5))
        text(expand.grid(x = 1:nrow(o3_distance), y = nrow(o3_distance):1), 
            labels = o3_distance, font = 2, cex = 0.6, col = "blue")
        par(xpd = TRUE)
        gradient.rect(xypos(3, 25)$x, xypos(3, 25)$y, xypos(7, 
            75)$x, xypos(7, 75)$y, nslices = 1001, col = colors, 
            gradient = "y")
        text(x = rep(xypos(7.5, 25)$x, 11), y = seq(xypos(7.5, 
            25)$y, xypos(7.5, 75)$y, length.out = 12), labels = c(as.character(round(seq(0, 
            1, length.out = 11), digits = 2)), "max"), pos = 4)
        text(x = xypos(2.3, 50)$x, y = xypos(2.3, 50)$y, labels = "Average Relative Difference", 
            srt = 90)
        mtext("Numbers in Ellipses Represent Distance Between Sites in km", 
            side = 1, line = 1)
        text(x = xypos(94, 78)$x, y = xypos(94, 78)$y, labels = expression(R^2), 
            font = 2)
        y <- 73
        trigger <- 0
        for (i in sqrt(c(1, seq(1, 0, by = -0.2)))) {
            if (trigger > 0) {
                plottype <- "l"
            }
            else {
                plottype <- "n"
            }
            subplot(plot(ellipse(i), axes = FALSE, type = plottype, 
                xlab = "", ylab = ""), xypos(90, y)$x, xypos(90, 
                y)$y, size = c(0.5, 0.5))
            text(x = xypos(94, y)$x, y = xypos(94, y)$y, labels = as.character(i^2), 
                font = 2)
            if (trigger > 0) {
                y <- y - 7
            }
            trigger <- 1
        }
        dev.off()
        setWinProgressBar(pb, value = 100, title = "Calculating", 
            label = "100% done")
        Sys.sleep(2)
        close(pb)
        return(1)
    }
    done <<- 0
    while (done != 3) {
        tclRequire("BWidget")
        dname <- ""
        dialog <- tktoplevel()
        tkwm.title(dialog, "Correlation Matrix Tool")
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
        yr <<- guiTextEntry(sframe = dialog, text = "Years to Include *", 
            default = "")
        tkgrid(yr$guiObject)
        tkgrid.configure(yr$guiObject, sticky = "nws")
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
                choosearea()
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
        flag <- 0
        while (flag == 0 & done != 3) {
            flag <- make.cormat()
            if (!flag) {
                tkmessageBox(title = "", message = "Chosen area has too few sites.\nPlease choose a larger area", 
                  icon = "error", type = "ok")
                choosearea()
            }
        }
        if (done != 3) {
            decide <- tktoplevel()
            tkwm.title(decide, "Decision Time!")
            make.more <- tkbutton(decide, text = "Make More\nPlots", 
                command = function() {
                  done <<- 0
                  tkdestroy(decide)
                })
            all.done <- tkbutton(decide, text = "All Done", command = function() {
                done <<- 3
                tkdestroy(decide)
            })
            tkgrid(tklabel(decide, text = "Your files are saved.\nDo you want to make more plots"))
            tkgrid(make.more, all.done)
            tkgrid.configure(all.done, sticky = "nw")
            tcl("focus", "-force", decide)
            tkwait.window(decide)
        }
    }
    message("Process completed")
}