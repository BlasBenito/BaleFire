#PLOT DIAGRAM
plotDiagram <- function(pollen.data, fire.data, fire.cluster, dendrogram, title, title.size, text.size, width, height, filename){
  
  #REQUIRED LIBRARIES
  require(ggplot2)
  require(viridis)
  
  #PREPARING DENDROGRAM
  if(dendrogram==TRUE){
  fire.cluster.ggplot <- hclust2ggplot(clust=fire.cluster,
                                       new.x=fire.data$age,
                                       min.y=0,
                                       max.y=1)
  }
  
  #GENERATING DATAFRAMES TO SHADE ALTERNATE CLUSTERING GROUPS
  group.stats = data.frame(xmin = tapply(fire.data$age, fire.data$clustering.zone, FUN = min), xmax =  tapply(fire.data$age, fire.data$clustering.zone, FUN = max), ymin = 0, ymax = Inf)
  even.groups = group.stats[which(!is.odd(1:nrow(group.stats))), ]
  odd.groups = group.stats[which(is.odd(1:nrow(group.stats))), ]
  
  #COMMON ELEMENTS
  #bit to add the shaded groups to the plots
  plot.even.groups <- geom_rect(data = even.groups, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray50", alpha = 0.2, inherit.aes = FALSE)
  plot.odd.groups <- geom_rect(data = odd.groups, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray80", alpha = 0.2, inherit.aes = FALSE)
  
  #theme to remove x axis completely
  no.x.axis <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())
  
  #common theme
  common.theme <- theme(text = element_text(size=text.size), legend.position="none", plot.margin = unit(c(0, 1, 0.1, 0), "cm"),  axis.text.y = element_text(size=text.size))
  
  #theme to remove y axis completely
  no.y.axis <- theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), text = element_text(size=text.size))
  
  #adjusting xlimits properly
  x.limits <- scale_x_continuous(expand = c(0, 0), limits = c(0,ceiling(max(pollen.data$age))), breaks=seq(0, ceiling(max(pollen.data$age)), by=1))
  
  
  #PLOTS
  
  #dendrogram "#3B528BFF"
  if(dendrogram==TRUE){
    
  p1 <- ggplot() +
    # plot.even.groups +
    # plot.odd.groups +
    geom_segment(data=fire.cluster.ggplot, aes(x=x, y=y, xend=xend, yend=yend), color="#3B528BFF") +
    theme(plot.title = element_text(size = title.size)) +
    ggtitle(title) +
    x.limits +
    no.x.axis +
    no.y.axis +
    common.theme +
    theme(axis.text.y = element_blank())
  
  p2 <- ggplot() +
    plot.even.groups +
    plot.odd.groups + 
    geom_segment(data=pollen.data, aes(xend=age, yend=0, x=age, y=ericaceae.par, color=ericaceae.par), size=1) + 
    scale_color_viridis(direction=-1) +
    ylab("Erica PAR \n (gr./cm2yr)") + 
    x.limits +
    no.x.axis +
    common.theme
  
  } else {
    
    #ericaceae par
    p2 <- ggplot() +
      plot.even.groups +
      plot.odd.groups + 
      geom_segment(data=pollen.data, aes(xend=age, yend=0, x=age, y=ericaceae.par, color=ericaceae.par), size=1) + 
      scale_color_viridis(direction=-1) +
      theme(plot.title = element_text(size = title.size)) +
      ggtitle(title) +
      ylab("Erica PAR \n (gr./cm2yr)") + 
      x.limits +
      no.x.axis +
      common.theme
  }
  
  #accumulation rate
  #   <- ggplot() +
  #   plot.even.groups +
  #   plot.odd.groups + 
  #   geom_segment(data=pollen.data, aes(xend=age, yend=0, x=age, y=acc.rate, color=acc.rate)) + 
  #   scale_colour_gradient(low = "gray90", high = "black") +
  #   ylab("Acc. rate \n (yr/cm)") + 
  #   ggtitle(title) +
  #   x.limits +
  #   no.x.axis +
  #   common.theme
  
   #charcoal
  p3 <- ggplot() +
    plot.even.groups +
    plot.odd.groups + 
    geom_segment(data=pollen.data, aes(xend=age, yend=0, x=age, y=charcoal.acc.rate, color=charcoal.acc.rate)) + 
    scale_color_viridis(direction=-1, option="inferno") +
    ylab("CHAR \n (part./cm2 yr)") + 
    x.limits +
    no.x.axis +
    common.theme
  
  
  #PLOTTING AND SAVING
  if(dendrogram==TRUE){
  print(cowplot::plot_grid(p1, p2, p3, align="v", ncol=1, rel_heights = c(1, rep(0.75, 2), 1.1)) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")))
  } else {
  print(cowplot::plot_grid(p2, p3, align="v", ncol=1, rel_heights = c(1, rep(0.75, 1), 1.1)) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")))
  }
  ggsave(filename = filename, width = width, height=height)
  
  
}
  


#FUNCTION TO SCALE DATA
scaleData=function(data, old.max, old.min, new.max, new.min){
  result=((data - old.min) / (old.max - old.min)) * (new.max - new.min) + new.min
  return(result)
}

#TO CHECK IF NUMBER IS ODD
is.odd <- function(x) x %% 2 != 0


#FUNCTION TO PREPARE CLUSTER AS DENDROGRAM TO BE PLOTTED IN GGPLOT2
hclust2ggplot <- function(clust, new.x, min.y, max.y){
  
  require(ggdendro)
  
  #clust to dendro
  clust.dendro <- dendro_data(clust)

  #extract segments to be plotted with ggplot2
  clust.segment <- segment(clust.dendro)
  
  #generating the reference data
  ref.data <- data.frame(old.x=1:max(clust.segment$xend), new.x)
  
  #joining together old x values
  old.x.values <- c(clust.segment$x, clust.segment$xend)
  
  #iterating through rows in ref.data to scale x values in dendrogram
  for(i in 2:nrow(ref.data)){
    
    #computing old and new min and max
    old.min <- ref.data[i-1, "old.x"]
    old.max <- ref.data[i, "old.x"]
    new.min <- ref.data[i-1, "new.x"]
    new.max <- ref.data[i, "new.x"]

    #selecting cases to rescale
    old.x.values.subset <- sort(unique(old.x.values[old.x.values >= old.min & old.x.values < old.max]))
    
    if(i==nrow(ref.data)){
      old.x.values.subset <- sort(unique(old.x.values[old.x.values >= old.min & old.x.values <= old.max])) 
    }
    
    #scaling cases
    new.values.subset <- scaleData(old.x.values.subset, old.max=old.max, old.min=old.min, new.max=new.max, new.min=new.min)
    
    #replacing values in data
    for(j in 1:length(old.x.values.subset)){
      clust.segment[clust.segment$x==old.x.values.subset[j], "x"]<-new.values.subset[j]
      clust.segment[clust.segment$xend==old.x.values.subset[j], "xend"]<-new.values.subset[j]
    }
    
  }#end of rescaling x values

  clust.segment$y <- scaleData(clust.segment$y, old.min=min(clust.segment$y), old.max=max(clust.segment$y), new.max=max.y, new.min=min.y)
  clust.segment$yend <- scaleData(clust.segment$yend, old.min=min(clust.segment$yend), old.max=max(clust.segment$yend), new.max=max.y, new.min=min.y)
  
  return(clust.segment)
}


#FUNCTION TO MERGE DATASETS AND INTERPOLATE THEM TO A REGULAR GRID
interpolateDatasets<-function(datasets.list, age.column.name, interpolation.time.step){
  
  #computing age ranges
  age.ranges<-sapply(datasets.list, FUN=function(x) range(x[, age.column.name]))
  #min of maximum ages
  min.age<-round(max(age.ranges[1,]), 1)
  #max of minimum ages
  max.age<-round(min(age.ranges[2,]), 1)
  
  #subsetting dataframes in list
  datasets.list<-lapply(datasets.list, function(x) x[x[, age.column.name] >= min.age & x[, age.column.name] <= max.age, ])
  
  #reference data
  reference.age <- seq(min.age, max.age, by=interpolation.time.step)
  
  #looping through datasets to interpolate
  for (dataset.to.interpolate in names(datasets.list)){
    
    #getting the dataset
    temp <- datasets.list[[dataset.to.interpolate]]
    
    #removing age from the colnames list
    colnames.temp <- colnames(temp)
    colnames.temp <- colnames.temp[which(colnames.temp != age.column.name)]
    
    #empty dataset to store interpolation
    temp.interpolated <- data.frame(age=reference.age)
    
    #iterating through columns
    for (column.to.interpolate in colnames.temp){
      
      #do not interpolate non-numeric columns
      if (is.numeric(temp[, column.to.interpolate])==FALSE){
        temp.interpolated[, column.to.interpolate]<-temp[, column.to.interpolate]
        next
      }
      
      #interpolation
      interpolation.formula <- as.formula(paste(column.to.interpolate, "~", age.column.name, sep=" "))
      
      #iteration through span values untill R-squared equals 1
      span.values=seq(0.01, -0.000001, by=-0.000001)
      for(span in span.values){
        
        interpolation.function = loess(interpolation.formula, data=temp, span=span, control=loess.control(surface="direct"))
        
        #check fit
        if(cor(interpolation.function$fitted, temp[, column.to.interpolate]) >=  0.99999999){break}
        
      }
      
      print(paste("Correlation between observed and interpolated data = ", cor(interpolation.function$fitted, temp[, column.to.interpolate]), sep=""))
      
      interpolation.function <- loess(interpolation.formula, data=temp, span=0.02, control=loess.control(surface="direct"))
      interpolation.result <- predict(interpolation.function, newdata=reference.age, se=FALSE)
      
      #constraining the range of the interpolation result to the range of the reference data
      interpolation.range<-range(temp[, column.to.interpolate])
      interpolation.result[interpolation.result < interpolation.range[1]] <- interpolation.range[1]
      interpolation.result[interpolation.result > interpolation.range[2]] <- interpolation.range[2]
      
      #putting the interpolated data back in place
      temp.interpolated[, column.to.interpolate]<-interpolation.result
      
    }#end of iteration through columns
    
    #removing the age column
    temp.interpolated[, age.column.name]=NULL
    
    #putting the data back in the list
    datasets.list[[dataset.to.interpolate]] <- temp.interpolated
    
  }#end of iterations through datasets
  
  #same rows?
  nrow.datasets<-sapply(datasets.list, FUN=function(x) nrow(x))
  if (length(unique(nrow.datasets))==1){
    
    #remove age from all dataframes
    datasets.list<-lapply(datasets.list, function(x) { x[, age.column.name] <- NULL; x })
    
    #put dataframes together
    output.dataframe <- do.call("cbind", datasets.list) #changes names
    
  } else {
    stop("Resulting datasets don't have the same number of rows, there's something wrong with something.")
  }
  
  #add reference.age
  output.dataframe <- data.frame(age=reference.age, output.dataframe)
  
  return(output.dataframe)
  
}


#GENERATES LAGGED CHARCOAL DATA AFTER EVERY ERICA SAMPLE
forwardLags <- function(lags, reference.data, data.to.lag){
  
  #df to store the lagged data
  lag.data <- data.frame(ericaceae.par=double(), charcoal.acc.rate=double(), lag=integer())
  
  #iterates through erica lines
  for (erica.case in nrow(erica):1){
    
    #take a line of the erica dataframe and replicate it as many times as lags are
    erica.value <- rep(erica[erica.case, "ericaceae.par"], max(lags))
    
    #get the age of the replicated line as age.erica
    erica.case.age<-erica[erica.case, "age"]
    erica.case.age.plus.lags<-round((erica.case.age - (0.01 * max(lags))), 3)
    
    #if beyond maximum age
    if (erica.case.age.plus.lags < min(char.interpolated$age)){break}
    
    #get from char.interpolated the lines with age > age.erica && age <= age.erica + lags
    char.temp<-char.interpolated[which(char.interpolated$age < erica.case.age & char.interpolated$age >= erica.case.age.plus.lags), "charcoal.acc.rate"]
    
    #put the data together (NOTE: lags have to be added in reverse)
    char.temp<-data.frame(ericaceae.par=erica.value, charcoal.acc.rate=char.temp, lag=rev(lags))
    names(char.temp)<-names(lag.data)
    
    #put them in the final table
    lag.data<-rbind(lag.data, char.temp)
    
  }#end of iterations
  
  #remove stuff we don't need
  rm(char.temp, erica.case, erica.case.age, erica.case.age.plus.lags, erica.value)
  
  #order by lag
  lag.data<-lag.data[order(lag.data$lag),]
  
  #standardize data
  #get lags column
  lags.column<-lag.data$lag
  
  #standardize
  lag.data<-scale(lag.data[, c("ericaceae.par", "charcoal.acc.rate")])
  
  #add lag
  lag.data<-data.frame(lag.data, lag=lags.column)
  
  #lag as factor
  lag.data$lag<-as.factor(lag.data$lag)
  
  return(lag.data)
}



#GENERATES LAGGED CHARCOAL DATA BEFORE EVERY ERICA SAMPLE
backwardLags <- function(lags, reference.data, data.to.lag){
  
  #df to store the lagged data
  lag.data <- data.frame(ericaceae.par=double(), charcoal.acc.rate=double(), lag=integer())
  
  #iterates through erica lines
  for (erica.case in 1:nrow(erica)){
    
    #take a line of the erica dataframe and replicate it as many times as lags are
    erica.value <- rep(erica[erica.case, "ericaceae.par"], max(lags))
    
    #get the age of the replicated line
    erica.case.age<-erica[erica.case, "age"]
    erica.case.age.plus.lags<-round((erica.case.age + (0.01 * max(lags))), 3)
    
    #if beyond maximum age
    if (erica.case.age.plus.lags > max(char.interpolated$age)){break}
    
    #get from char.interpolated the lines with age > age.erica && age <= age.erica + lags
    char.temp<-char.interpolated[which(char.interpolated$age > erica.case.age & char.interpolated$age <= erica.case.age.plus.lags), "charcoal.acc.rate"]
    
    #put the data together
    char.temp<-data.frame(ericaceae.par=erica.value, charcoal.acc.rate=char.temp, lag=lags)
    
    #put them in the final table
    lag.data<-rbind(lag.data, char.temp)
    
  }#end of iterations
  
  #remove stuff we don't need
  rm(char.temp, erica.case, erica.case.age, erica.case.age.plus.lags, erica.value)
  
  #order by lag
  lag.data<-lag.data[order(lag.data$lag),]
  
  #standardize data
  #get lags column
  lags.column<-lag.data$lag
  
  #standardize
  lag.data<-scale(lag.data[, c("ericaceae.par", "charcoal.acc.rate")])
  
  #add lag
  lag.data<-data.frame(lag.data, lag=lags.column)
  
  #lag as factor
  lag.data$lag<-as.factor(lag.data$lag)
  
  return(lag.data)
}

#FUNCTION TO INSERT MINOR TICKS IN GGPLOT
insert_minor <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}


#PLOT LAGS
plotLags <- function(forward.lagged.data, backward.lagged.data, lags, filename){
  
  pdf(filename, height=6, width=12)
  for (lag in lags){
    
    #temp data
    temp.future<-forward.lagged.data[forward.lagged.data$lag==lag, ]
    temp.past<-backward.lagged.data[backward.lagged.data$lag==lag, ]
    
    plot.past <- ggplot(data=temp.past, aes(x=ericaceae.par, y=charcoal.acc.rate)) + geom_point(shape=21, fill="gray50", color="black", size=3, alpha=0.5) +
      xlab("Erica PAR (pollen grains/cm2yr)") + 
      ylab("CHAR (particles/cm2 yr)") +
      ggtitle(paste("Erica vs. CHAR, backward lag ", lag, sep="")) +
      theme(text=element_text(size=12), plot.title=element_text(size = 16))
    
    plot.future <- ggplot(data=temp.future, aes(x=ericaceae.par, y=charcoal.acc.rate)) + geom_point(shape=21, fill="gray50", color="black", size=3, alpha=0.5) +
      xlab("Erica PAR (pollen grains/cm2yr)") + 
      ylab("CHAR (particles/cm2 yr)") +
      ggtitle(paste("Erica vs. CHAR, forward lag ", lag, sep="")) +
      theme(text=element_text(size=12), plot.title=element_text(size = 16))
    
    print(cowplot::plot_grid(plot.past, plot.future, align="h", ncol=2))
  }
  dev.off()

}


#COMPUTES GLS MODELS BETWEEN A RESPONSE AND A VARIABLE IN A LAGGED DATASET. PROVIDES A DATAFRAME
modelLagData <- function(model.formula, lagged.data){
  
  lags<-as.numeric(sort(unique(lagged.data$lag)))
  model.formula<-as.formula(model.formula)
  response <- all.vars(model.formula)[1]
  
  #list to store results
  results.list <- list()
  
  #fitting a model per lag
  for (lag in lags){
    results.list[[lag]] <- gls(model.formula, data=lagged.data[lagged.data$lag==lag,])
  }
  
  #list to store pseudo R2
  results.list.R2<-list()
  
  #computing pseudo R2 per lag
  for (lag in lags){
    results.list.R2[[lag]] <- cor(lagged.data[lagged.data$lag==lag, response], predict(results.list[[lag]]))^2
  }
  
  #gathering coefficients
  results.list.coef <- lapply(results.list, function(x){summary(x)$tTable[2,1]})
  
  #gathering the standard error of the coefficients
  results.list.coef.se <- lapply(results.list, function(x){summary(x)$tTable[2,2]})
  
  #gathering p-value of coefficients
  results.list.pvalue <- lapply(results.list, function(x){summary(x)$tTable[2,4]})
  
  #to data frame
  output.df <- as.data.frame(do.call(rbind, results.list.coef))
  output.df <- data.frame(lag=as.numeric(rownames(output.df)), Coefficient=output.df$V1)
  output.df[, "p-value"] <- as.data.frame(do.call(rbind, results.list.pvalue))$V1
  output.df$R2 <- as.data.frame(do.call(rbind, results.list.R2))$V1
  output.df$lag<-output.df$lag*10
  
  #se of coefficients
  output.df.se <- as.data.frame(do.call(rbind, results.list.coef.se))
  output.df.se$lower <- output.df$Coefficient - output.df.se$V1
  output.df.se$upper <- output.df$Coefficient + output.df.se$V1
  output.df.se$V1 <- NULL
  
  #to long format for plotting
  output.df.long <- gather(output.df, variable, value, 2:ncol(output.df))
  
  #adding the errors
  output.df.long$lower <- c(output.df.se$lower, output.df.long[output.df.long$variable=="p-value", "value"], output.df.long[output.df.long$variable=="R2", "value"])
  output.df.long$upper <- c(output.df.se$upper, output.df.long[output.df.long$variable=="p-value", "value"], output.df.long[output.df.long$variable=="R2", "value"])
  
  return(output.df.long)
}



#COMPUTING NULL
#randomizing lag column of both datasets 999 times and computing model with randomized data
modelRandomLagData <- function(lagged.data, model.formula, iterations){
  
  #getting response column
  response <- all.vars(as.formula(model.formula))[1]
  
  #list to store results
  results <- list()
  
  #computing one model per iteration
  for(i in 1:iterations){
    
    #randomization of the response column
    lagged.data[, response] <- lagged.data[sample(1:nrow(lagged.data)), response]
    
    #computing model and storing results in list
    results[[i]] <- modelLagData(model.formula=model.formula, lagged.data=lagged.data)
    
  }#end of iterations
  
  #preparing the data
  #getting lags and variable names
  null.model.df <- data.frame(lag=results[[1]]$lag, variable=results[[1]]$variable, stringsAsFactors = FALSE)
  
  #getting the value column of each dataframe in each list
  null.model.values <- sapply(results, `[[`, 3)
  
  null.model.df$value <- apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.5)
  null.model.df$upper <- apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.95)
  null.model.df$lower <- apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.05)
  
  return(null.model.df)
  
}


#FUNCTION TO PLOT MODELING RESULTS
plotModelOutput <- function(forward.results, forward.results.random, backward.results, backward.results.random, filename, width, height, title.size, text.size){
  
  #axes limits
  max.lag <- max(c(forward.results$lag, backward.results$lag))
  max.coefficient <- round(max(c(forward.results[forward.results$variable=="Coefficient", "value"], backward.results[backward.results$variable=="Coefficient", "value"], forward.results.random[forward.results.random$variable=="Coefficient", "upper"], backward.results.random[backward.results.random$variable=="Coefficient", "upper"])) + 0.1, 1)
  min.coefficient <- round(min(c(forward.results[forward.results$variable=="Coefficient", "value"], backward.results[backward.results$variable=="Coefficient", "value"],  forward.results.random[forward.results.random$variable=="Coefficient", "lower"], backward.results.random[backward.results.random$variable=="Coefficient", "lower"])) - 0.1, 1)
  max.R2 <- round(max(c(forward.results[forward.results$variable=="R2", "value"], backward.results[backward.results$variable=="R2", "value"])), 1)
  
  #separating pvalues
  # forward.results.pvalue <- forward.results[forward.results$variable=="p-value", c("lag", "value")]
  # backward.results.pvalue <- backward.results[backward.results$variable=="p-value", c("lag", "value")]
  
  #reference color scale
  viridis.colors <- viridis(10, option="D")
  
  #BACKWARD PLOT
  backward.plot.coefficient <- ggplot(data=subset(backward.results, variable=="Coefficient"), aes(x=lag, y=value)) +
    geom_ribbon(data=subset(backward.results.random, variable=="Coefficient"), aes(ymin=lower, ymax=upper), alpha=0.3, fill=viridis.colors[10]) + 
    geom_line(data=subset(backward.results.random, variable=="Coefficient"), aes(x=lag, y=value), alpha=0.6, color=viridis.colors[10], size=1.5) +
    geom_hline(yintercept=0, color="black", linetype=3) + 
    geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.5, fill=viridis.colors[5]) + 
    geom_line(size=1.5, color=viridis.colors[1]) + 
    ggtitle(expression("CHAR" %->% "Erica PAR")) + 
    theme(legend.position="none") + 
    xlab("") + 
    ylab("Standardized coefficient") + 
    scale_y_continuous(limits=c(min.coefficient, max.coefficient), breaks=seq(min.coefficient, max.coefficient, by=0.2)) +
    scale_x_reverse(limits=c(max.lag, 10), breaks=c(10, seq(100, max.lag, by=100))) +
    theme(axis.text = element_text(size=text.size), plot.title = element_text(size = title.size), plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x = element_blank())
  
  forward.plot.coefficient <- ggplot(data=subset(forward.results, variable=="Coefficient"), aes(x=lag, y=value, group=variable)) +
    geom_ribbon(data=subset(forward.results.random, variable=="Coefficient"), aes(ymin=lower, ymax=upper), alpha=0.3, fill=viridis.colors[10]) + 
    geom_line(data=subset(forward.results.random, variable=="Coefficient"), aes(x=lag, y=value), alpha=0.6, color=viridis.colors[10], size=1.5) +
    geom_hline(yintercept=0, color="black", linetype=3) + 
    geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.5, fill=viridis.colors[5]) + 
    geom_line(size=1.5, color=viridis.colors[1]) + 
    ggtitle(expression("Erica PAR" %->% "CHAR")) + 
    theme(legend.position="none") + 
    xlab("") + 
    ylab("") + 
    scale_y_continuous(limits=c(min.coefficient, max.coefficient), breaks=seq(min.coefficient, max.coefficient, by=0.2)) +
    scale_x_continuous(limits=c(10, max.lag), breaks=c(10, seq(100, max.lag, by=100))) +
    theme(axis.text = element_text(size=text.size), plot.title = element_text(size = title.size), plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.line.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank())
  
  backward.plot.R2 <- ggplot(data=subset(backward.results, variable=="R2"), aes(x=lag, y=value, group=variable)) +
    geom_ribbon(data=subset(backward.results.random, variable=="R2"), aes(ymin=lower, ymax=upper), alpha=0.3, fill=viridis.colors[10]) + 
    geom_line(data=subset(backward.results.random, variable=="R2"), aes(x=lag, y=value), alpha=0.6, color=viridis.colors[10], size=1.5) +
    geom_line(size=1.5, color=viridis.colors[1]) + 
    theme(legend.position="none") + 
    xlab("Years (before Erica samples)") + 
    ylab("Pseudo R squared") + 
    scale_y_continuous(limits=c(0, max.R2), breaks=seq(0, max.R2, by=0.1)) +
    scale_x_reverse(limits=c(max.lag, 10), breaks=c(10, seq(100, max.lag, by=100))) +
    theme(axis.text = element_text(size=text.size), axis.text.x = element_text(size=text.size), plot.title = element_text(size = title.size), plot.margin = unit(c(0.2, 0.5, 0, 0), "cm"))
  
  forward.plot.R2 <- ggplot(data=subset(forward.results, variable=="R2"), aes(x=lag, y=value, group=variable)) +
    geom_ribbon(data=subset(forward.results.random, variable=="R2"), aes(ymin=lower, ymax=upper), alpha=0.3, fill=viridis.colors[10]) + 
    geom_line(data=subset(forward.results.random, variable=="R2"), aes(x=lag, y=value), alpha=0.6, color=viridis.colors[10], size=1.5) +
    geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.5, fill=viridis.colors[5]) + 
    geom_line(size=1.5, color=viridis.colors[1]) + 
    theme(legend.position="none") + 
    xlab("Years (after Erica samples)") + 
    ylab("") + 
    scale_y_continuous(limits=c(0, max.R2), breaks=seq(0, max.R2, by=0.1)) +
    scale_x_continuous(limits=c(10, max.lag), breaks=c(10, seq(100, max.lag, by=100))) +
    theme(axis.text = element_text(size=text.size), axis.text.x = element_text(size=text.size), plot.title = element_text(size = title.size), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())
  
  first_col = plot_grid(backward.plot.coefficient, backward.plot.R2, ncol = 1, rel_heights = c(1, 1), align="v") + theme(plot.margin = unit(c(0.5, -1, 0.5, 0.5), "cm"))
  
  
  second_col = plot_grid(forward.plot.coefficient, forward.plot.R2, ncol = 1, rel_heights = c(1, 1), align="v") + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"))
  
  
  complete.plot = plot_grid(first_col, NULL, second_col, ncol = 3, rel_widths = c(1, 0, 1), align="h")
  
  print(complete.plot)
  
  ggsave(filename = filename, width=width, height=height)
  
}




