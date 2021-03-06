# load('F:/01 MetIMMS/00 data processing/190415 CCS consistency evaluation and outliner remove/test/raw_collected_V3.2_pos_metccs_test')
# metccs_fitting <- FitPowerFunction(raw_data = raw_collected_V3.2_pos_metccs,
#                                    id='id',
#                                    x = 'mz',
#                                    y = "X2016_ac_metccs")
#
# cred_para_set <- CredParaSetCalcu(fit_model = metccs_fitting$model,
#                                   fit_result = metccs_fitting$result)
#
# PiCalcu(x=metccs_fitting$result$x,
#         cred_para_set = cred_para_set)

#' @title FitPowerFunction
#' @description Power fitting
#' @author Zhiwei Zhou
#' @param raw_data data.frame, it need 3 columns at least
#' @param id character. column name of id
#' @param x character. column name of x
#' @param y character. column name of y
#' @export

setGeneric(name = 'FitPowerFunction',
           def = function(
             raw_data,
             id='id',
             x='x',
             y='y'
           ){
             temp_data <- raw_data

             colnames(temp_data)[which(colnames(temp_data)==x)] <- "x"
             colnames(temp_data)[which(colnames(temp_data)==y)] <- "y"
             colnames(temp_data)[which(colnames(temp_data)==id)] <- "id"

             fit_model <- nls(y ~ a*x^b,
                              data = temp_data,
                              start = list(a=1, b= 0.1))

             fit_result <- predict(fit_model, temp_data$x)

             fit_result <- data.frame(id=temp_data$id,
                                      x=temp_data$x,
                                      y=temp_data$y,
                                      pred_y=fit_result,
                                      stringsAsFactors = F)

             result <- list(model=fit_model,
                            result=fit_result)

             return(result)
           }
)

#' @title CredParaSetCalcu
#' @author Zhiwei Zhou
#' @description Calculate parameter set for CI and PI calucaltion
#' @param fit_model Power fitting model (Generated by FitPowerFunction)
#' @param fit_result Power fitting result (Generated by FitPowerFunction)
#' @export

setGeneric(name = 'CredParaSetCalcu',
           def = function(
             fit_model,
             fit_result
           ){
             temp <- summary(fit_model)
             a <- temp$parameters[1,1]
             b <- temp$parameters[2,1]
             n <- nrow(fit_result)

             mean_x <- mean(fit_result$x^b)
             mean_y <- mean(fit_result$y)
             sum_error_x_square <- sum(((fit_result$x)^b-mean_x)^2)
             sum_error_y_square <- sum((fit_result$y-mean_y)^2)
             sigma_square <- (sum_error_y_square - a*sum(((fit_result$x)^b-mean_x)*fit_result$y))/(n-2)
             sigma <- sigma_square^0.5

             result <- list(a=a,
                            b=b,
                            n=n,
                            mean_x=mean_x,
                            mean_y=mean_y,
                            sum_error_x_square=sum_error_x_square,
                            sum_error_y_square=sum_error_y_square,
                            sigma=sigma)

             return(result)
           }
)

#' @title CiCalcu
#' @author Zhiwei Zhou
#' @description Calculate confidence interval (CI) with input of x
#' @param x Power fitting model (Generated by FitPowerFunction)
#' @param cred_para_set Power fitting result (Generated by FitPowerFunction)
#' @param q vector of quantiles; Default: 0.995, representing 99% cridential
#' @export

setGeneric(name = 'CiCalcu',
           def = function(
             x,
             cred_para_set,
             q=0.995
           ){
             a <- cred_para_set[[1]]
             b <- cred_para_set[[2]]
             n <- cred_para_set[[3]]
             mean_x <- cred_para_set[[4]]
             mean_y <- cred_para_set[[5]]
             sum_error_x_square <- cred_para_set[[6]]
             sum_error_y_square <- cred_para_set[[7]]
             sigma <- cred_para_set[[8]]

             z_value <- qnorm(q, mean = 0, sd = 1)
             result <- sigma*(1/n+(x^b-mean_x)^2/sum_error_x_square)^0.5*z_value

             return(result)
           }
)


#' @title PiCalcu
#' @author Zhiwei Zhou
#' @description Calculate predictive interval (PI) with input of x
#' @param x Power fitting model (Generated by FitPowerFunction)
#' @param cred_para_set Power fitting result (Generated by FitPowerFunction)
#' @param q vector of quantiles; Default: 0.995, representing 99% cridential
#' @export

setGeneric(name = 'PiCalcu',
           def = function(
             x,
             cred_para_set,
             q=0.995
           ){
             a <- cred_para_set[[1]]
             b <- cred_para_set[[2]]
             n <- cred_para_set[[3]]
             mean_x <- cred_para_set[[4]]
             mean_y <- cred_para_set[[5]]
             sum_error_x_square <- cred_para_set[[6]]
             sum_error_y_square <- cred_para_set[[7]]
             sigma <- cred_para_set[[8]]

             z_value <- qnorm(q, mean = 0, sd = 1)
             result <- sigma*(1+1/n+(x^b-mean_x)^2/sum_error_x_square)^0.5*z_value

             return(result)
           }
)
