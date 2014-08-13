area.served <- function() 
{

    pb <- winProgressBar(title = "Area Served Tool Starting Up", 
        min = 0, max = 1, width = 500, label = "Loading data. Please wait")
    for (i in j <- c("census_00_centroids", "census_08")) {
        data(list = i, package = "netassess", envir = environment())
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
    data.test <- function(checkvar) {
        done <<- 0
        if (checkvar == "main3") {
            if (tclvalue(fname$object) == "" | dname == "" | 
                tclvalue(exist.sites$object) == "") {
                tt <- tkmessageBox(title = "'Doh!'", message = "One or more of the items marked by '*' is missing!", 
                  icon = "info", type = "ok")
            }
            else {
                done <<- 1
            }
        }
    }
    make.voronoi <- function() {
        pb <- winProgressBar("Making Voronoi Polygons", label = "0%", 
            width = 500, max = 100)
        if (tclvalue(new.sites$object) == "") {
            new <- NULL
        }
        else {
            new <- read.csv(tclvalue(new.sites$object), stringsAsFactors = FALSE, 
                colClasses = "character")
            new$type <- "new"
        }
        exist <- read.csv(tclvalue(exist.sites$object), stringsAsFactors = FALSE, 
            colClasses = "character")
        exist$type <- "existing"
        setWinProgressBar(pb, value = 10, title = "Making Voronoi Polygons", 
            label = "10% done")
        all.sites <- rbind(exist, new)
        names(all.sites) <- c("siteid", "longitude", "latitude", 
            "type")
        all.sites$longitude <- as.numeric(all.sites$longitude)
        all.sites$latitude <- as.numeric(all.sites$latitude)
        setWinProgressBar(pb, value = 20, title = "Making Voronoi Polygons", 
            label = "20% done")
        all.sites$coordlab <- paste(round(all.sites$longitude, 
            5), round(all.sites$latitude, 5), sep = ":")
        byvars <- split(all.sites, all.sites[, "coordlab"], drop = TRUE)
        all.sites <- lapply(byvars, function(x) {
            temp <- unique(x[, c("longitude", "latitude", "coordlab", 
                "type")])
            temp$siteid <- list(x$siteid)
            temp$coordlab <- NULL
            return(temp)
        })
        setWinProgressBar(pb, value = 30, title = "Making Voronoi Polygons", 
            label = "30% done")
        all.sites <- do.call("rbind", all.sites)
        all.sites <- cbind(all.sites, tolambert(all.sites[, c("longitude", 
            "latitude")]))
        vp <- tile.list(deldir(all.sites$lamx, all.sites$lamy))
        setWinProgressBar(pb, value = 40, title = "Making Voronoi Polygons", 
            label = "40% done")
        vp2 <- lapply(vp, function(z) {
            Polygon(as.matrix(cbind(c(z$x, z$x[1]), c(z$y, z$y[1]))))
        })
        names(vp2) <- do.call("rbind", lapply(all.sites$siteid, 
            paste, collapse = ":"))
        all.sites$ID <- names(vp2)
        setWinProgressBar(pb, value = 50, title = "Making Voronoi Polygons", 
            label = "50% done")
        vp <- lapply(names(vp2), function(z) {
            Polygons(list(vp2[[z]]), ID = z)
        })
        vp2 <- SpatialPolygons(vp, proj4string = CRS("+proj=lcc +lon_0=97w +lat_0=40n +lat_1=33n +lat_2=45n"))
        tract.in.vp <- overlay(census_00_centroids, vp2)
        setWinProgressBar(pb, value = 60, title = "Making Voronoi Polygons", 
            label = "60% done")
        census_fips <- as.character(census_00_centroids@data$FIPS)
        voronoi_ids <- sapply(slot(vp2, "polygons"), function(z) {
            slot(z, "ID")
        })
        setWinProgressBar(pb, value = 70, title = "Making Voronoi Polygons", 
            label = "70% done")
        voronoi <- data.frame(census_fips, voronoi_ids = voronoi_ids[tract.in.vp], 
            stringsAsFactors = F)
        voronoi <- merge(voronoi, census_00_centroids@data[, 
            c("FIPS", "POP2000")], by.x = "census_fips", by.y = "FIPS")
        voronoi <- merge(voronoi, census_08[, c("fips", "pop", 
            "SQMI")], by.x = "census_fips", by.y = "fips")
        setWinProgressBar(pb, value = 80, title = "Making Voronoi Polygons", 
            label = "80% done")
        voronoi_sum <- by(voronoi, list(voronoi_ids = voronoi$voronoi_ids), 
            function(z) {
                temp <- data.frame(voronoi_id = unique(z$voronoi_ids), 
                  tpop2000 = sum(z$POP2000), tpop = sum(z$pop), 
                  tarea = sum(z$SQMI), popden2000 = sum(z$POP2000)/sum(z$SQMI), 
                  popden = sum(z$pop)/sum(z$SQMI), stringsAsFactors = F)
                temp$tractfips <- list(as.character(z$census_fips))
                return(temp)
            })
        voronoi_sum <- do.call("rbind", voronoi_sum)
        setWinProgressBar(pb, value = 90, title = "Making Voronoi Polygons", 
            label = "90% done")
        voronoi_sum <- merge(voronoi_sum, data.frame(voronoi_id = voronoi_ids, 
            stringsAsFactors = FALSE), by = "voronoi_id", all.y = T)
        voronoi_sum <- merge(voronoi_sum, all.sites[, c("ID", 
            "longitude", "latitude", "type")], by.x = "voronoi_id", 
            by.y = "ID", all.x = TRUE)
        rownames(voronoi_sum) <- voronoi_sum$voronoi_id
        voronoi <- SpatialPolygonsDataFrame(vp2, voronoi_sum)
        voronoi <<- spTransform(voronoi, CRS("+proj=longlat +ellps=clrk66"))
        setWinProgressBar(pb, value = 100, title = "Making Voronoi Polygons", 
            label = "100% done")
        Sys.sleep(2)
        close(pb)
    }
    savefiles <- function() {
        toshp <- voronoi
        toshp$tractfips <- sapply(toshp$tractfips, paste, collapse = ",")
        writePolyShape(toshp, paste(dname, "/", tclvalue(fname$object), 
            sep = ""))
        write.csv(toshp@data, paste(dname, "/", tclvalue(fname$object), 
            ".csv", sep = ""))
        kml <- file(paste(dname, "/", tclvalue(fname$object), 
            ".kml", sep = ""), open = "w")
        cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", "<kml xmlns=\"http://earth.google.com/kml/2.1\">", 
            "<Document>", paste("\t<name>", tclvalue(fname$object), 
                "</name>", sep = ""), "\t<visibility>0</visibility>", 
            "\t<open>0</open>", "\t<Style id=\"vpolys-new\">", 
            "\t\t<PolyStyle>", "\t\t\t<color>00000000</color>", 
            "\t\t\t<colorMode>normal</colorMode>", "\t\t\t<fill>1</fill>", 
            "\t\t\t<outline>1</outline>", "\t\t</PolyStyle>", 
            "\t\t<LineStyle>", "\t\t\t<color>FF00FFFF</color>", 
            "\t\t\t<width>5</width>", "\t\t</LineStyle>", "\t\t<IconStyle>", 
            "\t\t\t<color>FFFF0000</color>", "\t\t\t<scale>1</scale>", 
            "\t\t\t<Icon>", "\t\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/star.png</href>", 
            "\t\t\t</Icon>", "\t\t</IconStyle>", "\t\t<LabelStyle>", 
            "\t\t\t<color>FF00FFFF</color>", "\t\t</LabelStyle>", 
            "\t</Style>", "\t<Style id=\"vpolys-exist\">", "\t\t<PolyStyle>", 
            "\t\t\t<color>00000000</color>", "\t\t\t<colorMode>normal</colorMode>", 
            "\t\t\t<fill>1</fill>", "\t\t\t<outline>1</outline>", 
            "\t\t</PolyStyle>", "\t\t<LineStyle>", "\t\t\t<color>FF00FFFF</color>", 
            "\t\t\t<width>5</width>", "\t\t</LineStyle>", "\t\t<IconStyle>", 
            "\t\t\t<color>FF0000FF</color>", "\t\t\t<scale>1</scale>", 
            "\t\t\t<Icon>", "\t\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/star.png</href>", 
            "\t\t\t</Icon>", "\t\t</IconStyle>", "\t\t<LabelStyle>", 
            "\t\t\t<color>FF00FFFF</color>", "\t\t</LabelStyle>", 
            "\t</Style>", sep = "\n", file = kml)
        garbage <- lapply(slot(voronoi, "polygons"), function(z) {
            cat("\t<Placemark>", paste("\t\t<name>", z@ID, "</name>", 
                sep = ""), "\t\t<description>", paste("\t\t\t<![CDATA[<hr /><table width=\"300\"><tr><td>Area: ", 
                round(voronoi@data[voronoi@data$voronoi_id == 
                  z@ID, "tarea"], digits = 0), " sqmi<br>2000 Population: ", 
                voronoi@data[voronoi@data$voronoi_id == z@ID, 
                  "tpop2000"], " (", round(voronoi@data[voronoi@data$voronoi_id == 
                  z@ID, "popden2000"], digits = 0), " people/sqmi)<br> 2008 Population: ", 
                voronoi@data[voronoi@data$voronoi_id == z@ID, 
                  "tpop"], " (", round(voronoi@data[voronoi@data$voronoi_id == 
                  z@ID, "popden"], digits = 0), " people/sqmi)<br><br>Census tract centroids within polygon:<br>", 
                sapply(voronoi@data[voronoi@data$voronoi_id == 
                  z@ID, "tractfips"], paste, collapse = "<br>"), 
                "</td></tr></table>]]>", sep = ""), "\t\t</description>", 
                sep = "\n", file = kml)
            if (voronoi@data[voronoi@data$voronoi_id == z@ID, 
                "type"] == "existing") {
                cat("\t\t<styleUrl>#vpolys-exist</styleUrl>", 
                  sep = "\n", file = kml)
            }
            else {
                cat("\t\t<styleUrl>#vpolys-new</styleUrl>", sep = "\n", 
                  file = kml)
            }
            cat("\t\t<MultiGeometry>", "\t\t<Polygon>", "\t\t\t<extrude>0</extrude>", 
                "\t\t\t<tessellate>0</tessellate>", "\t\t\t<altitudeMode>clampToGround</altitudeMode>", 
                "\t\t\t<outerBoundaryIs>", "\t\t\t<LinearRing>", 
                "\t\t\t<coordinates>", apply(slot(z, "Polygons")[[1]]@coords, 
                  1, function(xx) {
                    paste(xx[1], ",", xx[2], ",0", sep = "")
                  }), "\t\t\t</coordinates>", "\t\t\t</LinearRing>", 
                "\t\t\t</outerBoundaryIs>", "\t\t</Polygon>", 
                "\t\t<Point>", paste("\t\t\t<coordinates>", voronoi@data[voronoi@data$voronoi_id == 
                  z@ID, "longitude"], ",", voronoi@data[voronoi@data$voronoi_id == 
                  z@ID, "latitude"], ",0</coordinates>", sep = ""), 
                "\t\t</Point>", "\t\t</MultiGeometry>", "\t</Placemark>", 
                sep = "\n", file = kml)
        })
        cat("</Document>", "</kml>", sep = "\n", file = kml)
        close(con = kml)
        shell(paste("zip \"", gsub("/", "\\", dname, fixed = TRUE), 
            "\\", tclvalue(fname$object), ".kmz\" \"", gsub("/", 
                "\\", dname, fixed = TRUE), "\\", tclvalue(fname$object), 
            ".kml\"", sep = ""))
        shell(paste("del /q \"", gsub("/", "\\", dname, fixed = TRUE), 
            "\\", tclvalue(fname$object), ".kml\"", sep = ""))
        tkmessageBox(title = "", message = "Files saved.\nAll Done!", 
            icon = "info", type = "ok")
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
    main.win <- tktoplevel()
    tkwm.title(main.win, "Site Population Tool")
    exist.sites <- newguiFilename(main.win, text = "Existing Sites File *", 
        filter = "{{Comma Separated Values file} {.csv}}")
    new.sites <- newguiFilename(main.win, text = "New Sites File", 
        filter = "{{Comma Separated Values file} {.csv}}", )
    tkgrid(exist.sites$guiObject)
    tkgrid.configure(exist.sites$guiObject, sticky = "news")
    tkgrid(new.sites$guiObject)
    tkgrid.configure(new.sites$guiObject, sticky = "news")
    dframe <- guiFrame(sframe = main.win)
    tkgrid.configure(dframe, sticky = "news")
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
    fname <<- guiTextEntry(sframe = main.win, text = "File name *", 
        default = "")
    tkgrid(fname$guiObject)
    tkgrid.configure(fname$guiObject, sticky = "nesw")
    dframe1 <- tkframe(main.win)
    save.but <- tkbutton(main.win, text = "Make Voronoi polygons\nand save kmz, shp", 
        command = function() {
            data.test("main3")
            if (done == 1) {
                make.voronoi()
                savefiles()
            }
        })
    tkgrid(save.but, dframe1)
    tkgrid.configure(dframe1, sticky = "news")
    dframe2 <- tkframe(dframe1)
    cancel.but <- tkbutton(dframe2, text = "Cancel", command = function() {
        tkdestroy(main.win)
        stop("cancelled by user")
    })
    done.but <- tkbutton(dframe2, text = "Done", command = function() {
        tkdestroy(main.win)
    })
    tkgrid.configure(dframe2, sticky = "nes")
    tkgrid(done.but, cancel.but)
    tkgrid.configure(cancel.but, sticky = "nes")
    tkgrid.configure(cancel.but, sticky = "nes")
    tcl("focus", "-force", main.win)
    tkwait.window(main.win)
    rm(list = ls(pos = 1), pos = 1)
    gc(verbose = FALSE)
}