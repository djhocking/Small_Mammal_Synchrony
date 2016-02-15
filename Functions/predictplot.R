predict.plot.data.frame <- function(x,given,given.lab,layout,partial,
                                    type=c("prob","logit","probit"),
                                    identify.pred=F,scol="red",...) {
  type <- match.arg(type)
  resp <- response.var(x)
  pred <- predictor.vars(x)
  if(!missing(given) && !is.factor(given)) {
    # make given a factor
    b <- as.numeric(quantile(given,seq(0,1,length=5)))
    given <- rcut(given,b)
  }
  if(missing(layout)) {
    len <- length(pred) + !missing(given)
    layout <- c(1,len)
    if(len > 3) {
      layout[2] <- ceiling(sqrt(len))
      layout[1] <- ceiling(len/layout[2])
    }
  }
  opar <- par(mfrow=layout, mar=c(4.5,4,0,0.1))
  on.exit(par(opar))
  if(!missing(partial)) {
    pres <- data.frame(residuals(partial,type="partial"))
  }
  for(i in pred) {
    if(!missing(partial)) {
      x[[resp]] <- pres[[make.names(i)]]
      if(is.null(x[[resp]])) stop(paste("partial of",i,"not found"))
    }
    k <- !is.na(x[,i])
    if(is.factor(x[[i]]) && !is.ordered(x[[i]])) {
      x[[i]] <- sort.levels(x[[i]],as.numeric(x[[resp]]))
    }
    if(!missing(given) && is.factor(x[[i]])) {
      plot.new()
      xlim <- length(levels(x[[i]]))
      plot.window(xlim=c(0.5,xlim+0.5),ylim=range(x[[resp]]))
      axis(1,1:xlim,labels=levels(x[[i]]))
      axis(2)
      box()
      title(xlab=i, ylab=resp)
      cat(paste("jittering",i,"\n"))
    } else {
      if(is.factor(x[[resp]])) {
        if(type=="prob") {
          if(is.factor(x[[i]])) {
            mosaicplot(table(x[[i]], x[[resp]]), xlab=i, ylab=resp)
          } else {
            plot(x[[i]], x[[resp]], xlab=i, ylab=resp, ...)
          }
        }
      } else {
        plot(x[[i]], x[[resp]], xlab=i, ylab=resp, ...)
      }
    }
    if(missing(given) && !is.na(scol)) {
      if(is.factor(x[[resp]])) {
        if(length(levels(x[[resp]]))==2 && !is.factor(x[[i]])) {
        if(type=="prob") {
          lines(loprob(x[k,i], x[k,resp]), col=scol)
        } else {
          xy <- loprob(x[k,i], x[k,resp])
          p <- xy$y-1
          p <- switch(type,logit=log(p/(1-p)),probit=qnorm(p))
          xy$y <- p+1.5
          plot(xy,col=scol,type="l",xlab=i,ylab=type)
          points(x[[i]],2*as.numeric(x[[resp]])-3)
        }
        }
      } else {
        lines(lowess(x[k,i], x[k,resp]),col=scol)
      }
    }
    if((identify.pred == T) || (i %in% identify.pred)) {
      identify(x[k,i],x[k,resp],labels=rownames(x)[k])
    }
    if(!missing(given)) {
      lev <- levels(given)
      for(g in 1:length(lev)) {
        color <- ((g-1) %% 6) + 1
        val <- lev[g]
        k <- (given == val)
        if(is.factor(x[[i]])) {
          jitter <- (runif(length(x[k,i]))-0.5)/5
          points(as.numeric(x[k,i])+jitter, x[k,resp], col=color, ...)
        } else {
          points(x[k,i], x[k,resp], col=color, ...)
        }
        if(is.factor(x[[resp]])) {
          lines(loprob(x[k,i], x[k,resp]),col=color)
        } else {
          lines(lowess(x[k,i], x[k,resp]),col=color)
          #abline(lm(x[k,resp]~x[k,i]),col=color)
        }
      }
    }
  }
  if(!missing(given)) {
    # legend
    plot.new()
    if(!missing(given.lab)) title(xlab=given.lab)
    y <- cumsum(strheight(lev)+0.02)
    for(i in 1:length(lev)) {
      color <- ((i-1) %% 6) + 1
      val <- lev[i]
      text(0.5,0.75-y[i],val,col=color,adj=0.5) 
    }
  }
}
predict.plot.lm <- function(object,data,partial=F,...) {
  if(!partial) {
    if(missing(data)) {
      res <- residual.frame(object)
    } else {
      res <- residual.frame(object,data)
    }
    if(F) {
      expr <- match.call(expand = F)
      expr$... <- NULL
      expr[[1]] <- as.name("residual.frame")
      res <- eval(expr, parent.frame())
    }
    cat("plotting residuals\n")
    predict.plot.data.frame(res,...)
  } else {
    if(missing(data)) data <- model.frame(object)
    cat("plotting partial residuals\n")
    predict.plot.data.frame(data,partial=object,...)
  }
}
predict.plot.formula <- function(formula,data=parent.frame(),...) {
  # formula has givens?
  rhs <- formula[[3]]
  if(is.call(rhs) && (deparse(rhs[[1]]) == "|")) {
    # remove givens from formula
    given <- deparse(rhs[[3]])
    formula[[3]] <- rhs[[2]]
    if(is.environment(data)) g <- get(given,env=data)
    else g <- data[[given]]
    if(is.null(g)) 
      stop(paste("variable \"",given,"\" not found",sep=""))
    return(predict.plot.formula(formula,data,
                                given=g,given.lab=given,...))
  }
  if(F) {
    expr <- match.call(expand = F)
    expr$... <- NULL
    expr$na.action <- na.omit
    expr[[1]] <- as.name("model.frame.default")
    x <- eval(expr, parent.frame())
  } else {
    # formula has its own environment 
    x <- model.frame.default(formula,data,na.action=na.omit)
  }
  predict.plot.data.frame(x,...)
}
predict.plot <- function(object, ...) UseMethod("predict.plot")

step.up <- function(object) {
  resp <- response.var(object)
  pred <- predictor.vars(object)
  scope <- terms(formula(paste(resp,"~",paste(pred,collapse="*"))))
  step(object,scope)
}