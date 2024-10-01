#################################################################################################### .
#   Program/Function Name: data.preprocess
#   Author: Ya Wang
#   Description: Preprocess input dataset to include required variables and output as a matrix
#   Change History:
#   Last Modified Date: 09/18/2024
#################################################################################################### .
#' @name data.preprocess
#' @title data.preprocess
#' @description { Preprocess input dataset to include required variables and output as a matrix }
#' @param data input dataset (a data frame)
#' @param resp the name of the response variable, for survival endpoint, it is the time-to-event variable
#' @param event (for survival endpoint only) the name of the event variable, which takes value 1 for event and 0 for censored, or NULL for other types of endpoint
#' @param trt the name of the treatment variable, which takes value 1 for experimental treatment and 0 for control treatment
#' @param stratcov a vector containing the name(s) of stratification factors, or NULL
#' @param basecov.cont a vector containing the name(s) of continuous baseline covariates, or NULL
#' @param basecov.cat a vector containing the name(s) of categorical baseline covariates, or NULL
#' @return A matrix with required variables in specific order

data.preprocess <- function( data, resp, event = NULL, trt, stratcov = NULL, basecov.cont = NULL, basecov.cat = NULL)
{

     dfData <- as.data.frame( data )

     dfData[ , trt ] <- as.factor( dfData[ , trt ] )
     dfData[ , trt ] <- droplevels( dfData[ , trt ], exclude = c( "", NA ) )
     dfData[ , trt ] <- as.numeric( dfData[ , trt ] ) - 1

    if( !is.null( stratcov ) )
    {
        if( length( stratcov ) > 1 )
        {
            dfData$stratum <- apply( dfData[ , stratcov ] , 1 , paste , collapse = "_" )
        }
        else
        {
            dfData$stratum <- dfData[ , stratcov ]
        }

        dfData$stratum <- as.factor( dfData$stratum )

        checkmissingstratum <- apply( dfData[ , stratcov, drop = F ], 1, function( i ){
            sum( is.na( i ) | i == "" )
        } )
        dfData$stratum[ which( checkmissingstratum > 0 ) ] <- NA

        dfData$stratum <- droplevels( dfData$stratum, exclude = c( "", NA ) )
        dfData$stratum <- as.numeric( dfData$stratum )

        varNames <- c( resp, event, trt, "stratum", basecov.cont )
    }
    else
    {
        varNames <- c( resp, event, trt, basecov.cont )
    }


    if( !is.null( basecov.cat ) )
    {
        for( i in 1:length( basecov.cat ) )
        {

            cov.name  <- basecov.cat[ i ]
            cov.class <- class( dfData[ , cov.name ] )

            if( cov.class != "factor" )
            {
                dfData[ , cov.name ] <- as.factor( dfData[ , cov.name ] )
            }


            cov.levels <- levels( dfData[ , cov.name ] )
            cov.levels <- cov.levels[ cov.levels!="" ]
            n.levels   <- length( cov.levels )

            if( n.levels > 1 ) # convert categorical covariate with more then 2 levels to dummy variables
            {
                for( j in 1:( n.levels - 1 ) )
                {
                    dfData[ , paste( cov.name, j, sep = "_" ) ] <- ifelse( is.na( dfData[ , cov.name ] ) | dfData[ , cov.name ] == "", NA, ifelse( dfData[ , cov.name ] == cov.levels[ j + 1 ], 1, 0 ) )

                    varNames <- c( varNames, paste( cov.name, j, sep = "_" ) )
                }
            }
        }
    }




    dfData     <- dfData[ , c( varNames ) ]
    matData    <- as.matrix( dfData )

    return( matData )

}

