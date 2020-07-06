corrPartClass <- R6::R6Class(
    "corrPartClass",
    inherit = corrPartBase,
    private = list(
        #### Init + run functions ----
        .init=function() {

            private$.initCorrTable()

        },
        .run=function() {

            results <- private$.compute()

            private$.populateCorrTable(results)

        },

        #### Compute results ----
        .compute = function() {

            vars <- self$options$vars
            nVars <- length(vars)
            hyp <- self$options$hypothesis
            controls <- self$options$controls
            type <- self$options$type

            results <- list()
            if (nVars > 1) {
                for (i in 1:nVars) {
                    rowVar <- vars[[i]]
                    for (j in 1:nVars) {
                        if (j >= i) { next }
                        colVar <- vars[[j]]
                        resLst <- private$.test(rowVar, colVar, controls, type, hyp)
                        results[[rowVar]][[colVar]] <- resLst[, 1]
                        if (type == 'semi') { results[[colVar]][[rowVar]] <- resLst[, 2] }
                    }
                }
            }
            return(results)
        },

        #### Init tables/plots functions ----
        .initCorrTable = function() {

            matrix <- self$results$matrix
            vars <- self$options$vars
            nVars <- length(vars)
            type <- self$options$type

            for (i in seq_along(vars)) {
                var <- vars[[i]]

                matrix$addColumn(name=paste0(var, '[r]'), title=var,
                                 type='number', format='zto', visible='(pearson)')
                matrix$addColumn(name=paste0(var, '[rp]'), title=var,
                                 type='number', format='zto,pvalue', visible='(pearson && sig)')
                matrix$addColumn(name=paste0(var, '[rho]'), title=var,
                                 type='number', format='zto', visible='(spearman)')
                matrix$addColumn(name=paste0(var, '[rhop]'), title=var,
                                 type='number', format='zto,pvalue', visible='(spearman && sig)')
                matrix$addColumn(name=paste0(var, '[tau]'), title=var,
                                 type='number', format='zto', visible='(kendall)')
                matrix$addColumn(name=paste0(var, '[taup]'), title=var,
                                 type='number', format='zto,pvalue', visible='(kendall && sig)')
                matrix$addColumn(name=paste0(var, '[n]'), title=var,
                                 type='integer', visible='(n)')
            }

            for (i in seq_along(vars)) {

                var <- vars[[i]]
                values <- list()

                if (type != 'semi') {
                    for (j in seq(i, nVars)) {
                        v <- vars[[j]]
                        values[[paste0(v, '[r]')]] <- ''
                        values[[paste0(v, '[rp]')]] <- ''
                        values[[paste0(v, '[rho]')]] <- ''
                        values[[paste0(v, '[rhop]')]] <- ''
                        values[[paste0(v, '[tau]')]] <- ''
                        values[[paste0(v, '[taup]')]] <- ''
                        values[[paste0(v, '[n]')]] <- ''
                    }
                }

                values[[paste0(var, '[r]')]] <- '\u2014'
                values[[paste0(var, '[rp]')]] <- '\u2014'
                values[[paste0(var, '[rho]')]] <- '\u2014'
                values[[paste0(var, '[rhop]')]] <- '\u2014'
                values[[paste0(var, '[tau]')]] <- '\u2014'
                values[[paste0(var, '[taup]')]] <- '\u2014'
                values[[paste0(var, '[n]')]] <- '\u2014'

                matrix$setRow(rowKey=var, values)
            }

            hyp <- self$options$hypothesis
            flag <- self$options$flag
            controls <- self$options$controls
            nControls <- length(controls)

            if (length(controls) > 0) {
                matrix$setNote('controls', jmvcore::format('controlling for {}', listItems(controls)))
            }

            if (hyp == 'pos') {
                matrix$setNote('hyp', 'H\u2090 is positive correlation')
                if (flag)
                    matrix$setNote('flag', '* p < .05, ** p < .01, *** p < .001, one-tailed')
            }
            else if (hyp == 'neg') {
                matrix$setNote('hyp', 'H\u2090 is negative correlation')
                if (flag)
                    matrix$setNote('flag', '* p < .05, ** p < .01, *** p < .001, one-tailed')
            }
            else {
                matrix$setNote('hyp', NULL)
                if (flag)
                    matrix$setNote('flag', '* p < .05, ** p < .01, *** p < .001')
            }

            if ( ! flag)
                matrix$setNote('flag', NULL)

            if (type == 'part' && nControls > 0) {
                titleMatrix <- "Partial Correlation"
            } else if (type == 'semi' && nControls > 0) {
                titleMatrix <- "Semipartial Correlation"
            } else {
                titleMatrix <- "Correlation"
            }

            pearson <- self$options$pearson
            spearman <- self$options$spearman
            kendall <- self$options$kendall
            n <- self$options$n
            sig <- self$options$sig

            if ( ! sum(pearson, spearman, kendall) > 1  && ! n && ! sig) {
                if (pearson)
                    titleMatrix <- jmvcore::format("{} - Pearson's r", titleMatrix)
                else if (spearman)
                    titleMatrix <- jmvcore::format("{} - Spearman's rho", titleMatrix)
                else if (kendall)
                    titleMatrix <- jmvcore::format("{} - Kendall's Tau B", titleMatrix)
            }

            matrix$setTitle(titleMatrix)
        },

        #### Populate tables functions ----
        .populateCorrTable = function(results) {

            matrix <- self$results$matrix
            vars <- self$options$vars
            nVars <- length(vars)
            flag <- self$options$flag
            type <- self$options$type

            if (nVars > 1) {
                for (i in 1:nVars) {
                    rowVarName <- vars[[i]]
                    for (j in 1:nVars) {

                        if (i == j || (type != 'semi' && j > i))
                            next

                        values <- list()
                        colVarName <- vars[[j]]
                        result = results[[rowVarName]][[colVarName]]

                        values[[paste0(colVarName, '[r]')]] <- result$r
                        values[[paste0(colVarName, '[rp]')]] <- result$rp
                        values[[paste0(colVarName, '[rho]')]] <- result$rho
                        values[[paste0(colVarName, '[rhop]')]] <- result$rhop
                        values[[paste0(colVarName, '[tau]')]] <- result$tau
                        values[[paste0(colVarName, '[taup]')]] <- result$taup

                        values[[paste0(colVarName, '[n]')]] <- result$n

                        matrix$setRow(rowNo=i, values)

                        if (flag) {
                            if (result$rp < .001)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[r]'), '***')
                            else if (result$rp < .01)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[r]'), '**')
                            else if (result$rp < .05)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[r]'), '*')

                            if (result$rhop < .001)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[rho]'), '***')
                            else if (result$rhop < .01)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[rho]'), '**')
                            else if (result$rhop < .05)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[rho]'), '*')

                            if ( ! self$options$kendall)
                            {}  # do nothing
                            else if (result$taup < .001)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[tau]'), '***')
                            else if (result$taup < .01)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[tau]'), '**')
                            else if (result$taup < .05)
                                matrix$addSymbol(rowNo=i, paste0(colVarName, '[tau]'), '*')
                        }
                    }
                }
            }
        },

        #### Helper functions ----
        .cleanData = function(var1, var2) {

            dataRaw <- self$data
            controls <- self$options$controls
            type <- self$options$type

            data <- list()
            data[[var1]] <- jmvcore::toNumeric(dataRaw[[var1]])
            data[[var2]] <- jmvcore::toNumeric(dataRaw[[var2]])

            for (control in controls) {
                data[[control]] <- jmvcore::toNumeric(dataRaw[[control]])
            }

            attr(data, 'row.names') <- seq_len(length(data[[1]]))
            attr(data, 'class') <- 'data.frame'
            data <- jmvcore::naOmit(data)

            return(data)
        },
        .test = function(var1, var2, controls, type, hyp) {

            data <- private$.cleanData(var1, var2)
            var1 <- data[[var1]]
            var2 <- data[[var2]]

            n <- length(data[[1]])
            d <- 1 + as.integer(type == 'semi')
            c <- length(controls)
            if (c == 0) { rdta <- cbind(var1, var2) } else { rdta = cbind(var1, var2, data[, controls]) }
            results <- replicate(d, list(r = NaN, rp = 1, rho = NaN, rhop = 1, tau = NaN, taup = 1, n = n))

            suppressWarnings({
                for (method in c('pearson', 'spearman', 'kendall')) {
                    # calculate 'kendall' only if the respective tick box is set
                    if (method == 'kendall' && self$options$kendall == F) { next }
                        
                    try({
                        # the following code took some inspiration from the ppcor-package (https://cran.r-project.org/web/packages/ppcor)
                        # calculate covariance matrix and invert it
                        cvx <- cov(rdta, method = method)
	                    if (det(cvx) < .Machine$double.eps) { icvx <- ginv(cvx) } else { icvx <- solve(cvx) }
	                
	                    # calculation of the partial correlation coefficient 
                        if (type == 'part') {
                            results[ifelse(method == 'pearson', 'r', ifelse(method == 'spearman', 'rho', 'tau')), 1] <- -cov2cor(icvx)[1, 2]
                        # calculation of the semi-partial correlation coefficients
                        # this is more complex, to speed up processing both row and column variables
                        # are calculated at once and given back in a list with two entries
                        } else {
                            spcor <- -cov2cor(icvx) / sqrt(diag(cvx)) / sqrt(abs(diag(icvx) - t(t(icvx ^ 2) / diag(icvx))))
                            results[ifelse(method == 'pearson', 'r', ifelse(method == 'spearman', 'rho', 'tau')), 1] <- spcor[1, 2]
                            results[ifelse(method == 'pearson', 'r', ifelse(method == 'spearman', 'rho', 'tau')), 2] <- spcor[2, 1]
                        }
                
	                # assign p-value 
                        # repeated depending on whether it is a semi-partial correlation (d = 2) or not (d = 1)
                        for (i in 1:d) {            
                            # calculate the p-values for undirected ('corr') and directed hypotheses ('pos', 'neg')
                            # if the correlation coefficient doesn't match the proposed direction, the value is not set and remains NA
                            if (hyp == 'corr') {
                                if      (method == 'pearson')  { results[, i]$rp   <- 2 *    pt(-abs(results[, i]$r   * sqrt((n - 2 - c) / (1 - results[, i]$r   ^ 2))), (n - 2 - c)) }
                                else if (method == 'spearman') { results[, i]$rhop <- 2 *    pt(-abs(results[, i]$rho * sqrt((n - 2 - c) / (1 - results[, i]$rho ^ 2))), (n - 2 - c)) }
                                else if (method == 'kendall')  { results[, i]$taup <- 2 * pnorm(-abs(results[, i]$tau / sqrt(2 * (2 * (n - c) + 5) / (9 * (n - c) * (n - 1 - c)))))   }
                            } else if ((hyp == 'neg' && results[, i]$r   < 0) || (hyp == 'pos' && results[, i]$r   > 0)) {
                                if      (method == 'pearson')  { results[, i]$rp   <-        pt(-abs(results[, i]$r   * sqrt((n - 2 - c) / (1 - results[, i]$r   ^ 2))), (n - 2 - c)) }
                                else if (method == 'spearman') { results[, i]$rhop <-        pt(-abs(results[, i]$rho * sqrt((n - 2 - c) / (1 - results[, i]$rho ^ 2))), (n - 2 - c)) }
                                else if (method == 'kendall')  { results[, i]$taup <-     pnorm(-abs(results[, i]$tau / sqrt(2 * (2 * (n - c) + 5) / (9 * (n - c) * (n - 1 - c)))))   }
                            }
                        } # for
                    }) # try
                } # for
            }) # suppressWarnings
            results
        })
)
