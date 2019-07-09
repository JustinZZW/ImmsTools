#' @title RunBioTransformer
#' @author Zhiwei Zhou
#' @description R interface to run BioTransformer in CMD. Please refer Djoumbou-Feunang Y, Fiamoncini J, de la Fuente AG, Manach C, Greiner R, and Wishart DS; BioTransformer: A Comprehensive Computational Tool for Small Molecule Metabolism Prediction and Metabolite Identification; Journal of Cheminformatics 201911:2; DOI: 10.1186/s13321-018-0324-5.
#' @param pgm_path the path of biotransformer program. Default: 'I:/software/Tools/BioTransformer-master'
#' @param smiles the smiles of structure
#' @param output_file the csv name of output
#' @param task 'pred', 'cid'. prediction or compound identification. Default: pred
#' @param bt_type 'ecbased', 'cyp450', 'phaseII', 'hgut', 'superbio', 'allHuman', 'envimicro'. Default: 'ecbased'. Ref: https://github.com/JustinZZW/BioTransformer
#' @param str_type Default: '-ismi'
#' @param output_format Default: '-ocsv'
#' @param n_step number for reaction. Default: 1
#' @export
#' @examples
#' RunBioTransformer(pgm_path='I:/software/Tools/BioTransformer-master',
#' smiles = "CC(C)C1=CC=C(C)C=C1O",
#' output_file='I:/software/Tools/BioTransformer-master/test.csv',
#' n_step=1)


setGeneric(name = 'RunBioTransformer',
           def = function(
             pgm_path,
             smiles,
             output_file,
             task=c('pred', 'cid'),
             bt_type=c('ecbased', 'cyp450', 'phaseII', 'hgut', 'superbio', 'allHuman', 'envimicro'),
             str_type=c('-ismi'),
             output_format=c('-ocsv'),
             n_step=1
           ){
             task <- match.arg(task)
             bt_type <- match.arg(bt_type)
             str_type <- match.arg(str_type)
             output_format <- match.arg(output_format)

             cmd1 <- paste0(stringr::str_split(pgm_path, pattern = '\\:/')[[1]][1], ':')
             cmd2 <- paste0('cd ', pgm_path)
             cmd3 <- paste('java -jar',
                           # pgm_path,
                           # file.path(pgm_path, 'biotransformer-1-0-8.jar'),
                           'biotransformer-1-0-8.jar',
                           '-k',
                           task,
                           '-b',
                           bt_type,
                           str_type,
                           smiles,
                           output_format,
                           output_file,
                           '-s',
                           n_step,
                           sep = ' ')

             # run cmd sequensly
             cmd <- paste(cmd1, cmd2, cmd3, sep = ' && ')
             shell(cmd)
           }
)

