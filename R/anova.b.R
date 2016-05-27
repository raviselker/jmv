
AnovaClass <- R6::R6Class(
    "AnovaClass",
    inherit=silkycore::Analysis,
    private=list(
        .model=NA,
        .run=function() {
            
            dependentName <- self$options$get('dependent')
            fixedFactors <- self$options$get('fixedFactors')
            
            if (is.null(dependentName) || is.null(fixedFactors))
                return()
            
            data <- self$data
            
            for (factorName in fixedFactors)
                data[[factorName]] <- as.factor(data[[factorName]])
            
            data[[dependentName]] <- as.numeric(data[[dependentName]])
            
            for (contrast in self$options$get('contrasts')) {
                base::options(contrasts=c("contr.sum","contr.poly"))
                levels <- base::levels(data[[contrast$var]])
                stats::contrasts(data[[contrast$var]]) <- private$.createContrasts(levels, contrast$type)
            }
            
            factors <- paste0('`', self$options$get('fixedFactors'), '`', collapse='*')
            formula <- paste0('`', self$options$get('dependent'), '`~', factors)
            
            formula <- stats::as.formula(formula)
            private$.model <- stats::aov(formula, data)
            results <- car::Anova(private$.model, type='II', singular.ok=TRUE)
            
            anovaTable <- self$results$get('anova')
            rowCount <- dim(results)[1]
            rowNames <- dimnames(results)[[1]]
            
            for (i in 1:rowCount) {
                rowName <- rowNames[i]
                
                ss <- results[i,'Sum Sq']
                df <- results[i,'Df']
                ms <- results[i,'Sum Sq'] / results[i,'Df']
                F  <- results[i,'F value']
                p  <- results[i,'Pr(>F)']
                
                if ( ! is.finite(ss))
                    ss <- NA
                if (df == 0)
                    df <- NA
                if ( ! is.finite(ms))
                    ms <- NA
                if ( ! is.finite(F))
                    F <- NA
                if ( ! is.finite(p))
                    p <- NA
                
                tableRow <- list(ss=ss, df=df, ms=ms, F=F, p=p)
                anovaTable$addRow(rowKey=rowName, tableRow)
            }
            
            private$.populateContrasts(data)
            private$.populateLevenes(data)
        },
        .populateContrasts=function(data) {
            
            contrResults <- stats::summary.lm(private$.model)[["coefficients"]]
            contrasts <- self$options$get('contrasts')
            
            for (contrast in contrasts) {
                
                var <- contrast$var
                type <- contrast$type
                levels <- base::levels(data[[var]])
                labels <- private$.contrastLabels(levels, type)
                
                table <- self$results$get('contrasts')$get(var)
                
                for (i in seq_along(labels)) {
                    label <- labels[[i]]
                    name <- paste0(var, i)
                    table$addRow(rowKey=i, list(
                        contrast=label,
                        est=contrResults[name, "Estimate"],
                        se=contrResults[name, "Std. Error"],
                        t=contrResults[name, "t value"],
                        p=contrResults[name, "Pr(>|t|)"]
                    ))
                }
            }
        },
        .populateLevenes=function(data) {
            
            if ( ! self$options$get('homoTests'))
                return()
            
            dep <- self$options$get('dependent')
            factors <- self$options$get('fixedFactors')
            rhs <- paste0(factors, collapse=':')
            formula <- as.formula(paste0(dep, '~', rhs))
            
            result <- car::leveneTest(formula, data, center="mean")
            
            table <- self$results$get('assump')$get('eqVar')
            
            print(result)
            
            table$setRow(rowNo=1, values=list(
                F=result[1,'F value'],
                df1=result[1,'Df'],
                df2=result[2,'Df'],
                p=result[1,'Pr(>F)']))
        },
        .contrastLabels=function(levels, type) {
            
            nLevels <- length(levels)
            labels <- list()
            
            if (length(levels) <= 1) {
                
                # do nothing
                
            } else if (type == 'simple') {
            
                for (i in seq_len(nLevels-1))
                    labels[[i]] <- paste(levels[i+1], '-', levels[1])
                
            } else if (type == 'deviation') {
                
                all <- paste(levels, collapse=', ')
                for (i in seq_len(nLevels-1))
                    labels[[i]] <- paste(levels[i+1], '-', all)
                
            } else if (type == 'difference') {
                
                for (i in seq_len(nLevels-1)) {
                    rhs <- paste0(levels[1:i], collapse=', ')
                    labels[[i]] <- paste(levels[i + 1], '-', rhs)
                }
                
            } else if (type == 'helmert') {
                
                for (i in seq_len(nLevels-1)) {
                    rhs <- paste(levels[(i+1):nLevels], collapse=', ')
                    labels[[i]] <- paste(levels[i], '-', rhs)
                }
                
            } else if (type == 'polynomial') {
                
                names <- c('linear', 'quadratic', 'cubic', 'quartic', 'quintic', 'sextic', 'septic', 'octic')
                
                for (i in seq_len(nLevels-1)) {
                    if (i <= length(names)) {
                        labels[[i]] <- names[i]
                    } else {
                        labels[[i]] <- paste('degree', i, 'polynomial')
                    }
                }
            }
            
            labels
        },
        .createContrasts=function(levels, type) {
            
            nLevels <- length(levels)
            
            if (type == 'simple') {
                
                dummy <- contr.treatment(levels)
                dimnames(dummy) <- NULL
                coding <- matrix(rep(1/nLevels, prod(dim(dummy))), ncol=nLevels-1)
                contrast <- (dummy - coding) * nLevels
                
            } else if (type == 'deviation') {
                
                contrast <- matrix(0, nrow=nLevels, ncol=nLevels-1)
                
                for (i in seq_len(nLevels-1)) {
                    contrast[i+1, i] <- -1
                    contrast[1, i] <- 1
                }
                
            } else if (type == 'difference') {
                
                contrast <- stats::contr.helmert(levels)
                dimnames(contrast) <- NULL
                
            } else if (type == 'helmert') {
                
                contrast <- matrix(0, nrow=nLevels, ncol=nLevels-1)
                
                for (i in seq_len(nLevels-1)) {
                    p <- (1 / (nLevels - i + 1))
                    contrast[i,i] <- p * (nLevels - i)
                    contrast[(i+1):nLevels,i] <- -p
                }
                
                contrast
                
            } else if (type == 'polynomial') {
                
                contrast <- stats::contr.poly(levels)
                dimnames(contrast) <- NULL
            }
            
            contrast
        })
)
