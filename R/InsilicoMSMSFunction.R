#' @title MetFragMatch
#' @author Zhiwei Zhou
#' @param  DatabaseSearchRelativeMassDeviation Default: 5ppm
#' @param FragmentPeakMatchAbsoluteMassDeviation Default: 0.001Da
#' @param FragmentPeakMatchRelativeMassDeviation Default: 5.0ppm
#' @param MetFragDatabaseType 'PubChem'
#' @param NeutralPrecursorMass 253.966126
#' @param PeakList matrix with 2 columns, 1st: mz, 2nd: intensity
#' @return a data.frame with scores and candidate properties
#' @export
#' @example
#' test_spec <- matrix(c(
#' 90.97445, 681,
#' 106.94476, 274,
#' 110.02750, 110,
#' 115.98965, 95,
#' 117.98540, 384,
#' 124.93547, 613,
#' 124.99015, 146,
#' 125.99793, 207,
#' 133.95592, 777,
#' 143.98846, 478,
#' 144.99625, 352,
#' 146.00410, 999,
#' 151.94641, 962,
#' 160.96668, 387,
#' 163.00682, 782,
#' 172.99055, 17,
#' 178.95724, 678,
#' 178.97725, 391,
#' 180.97293, 999,
#' 196.96778, 720,
#' 208.96780, 999,
#' 236.96245, 999,
#' 254.97312, 999), ncol=2, byrow=TRUE)
#' MetFragMatch(DatabaseSearchRelativeMassDeviation = 5,
#'              FragmentPeakMatchAbsoluteMassDeviation = 0.001,
#'              FragmentPeakMatchRelativeMassDeviation = 5,
#'              MetFragDatabaseType = 'PubChem',
#'              NeutralPrecursorMass = 253.966126,
#'              PeakList = test_spec)


setGeneric(name = 'MetFragMatch',
           def = function(
             DatabaseSearchRelativeMassDeviation = 5, # ppm
             FragmentPeakMatchAbsoluteMassDeviation = 0.001, # Da
             FragmentPeakMatchRelativeMassDeviation = 5.0, # ppm
             MetFragDatabaseType = 'PubChem',
             NeutralPrecursorMass = 253.966126,
             PeakList
           ){

             # creat parameter set
             settingsObject<-list()
             settingsObject[["DatabaseSearchRelativeMassDeviation"]] <- DatabaseSearchRelativeMassDeviation
             settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]] <- FragmentPeakMatchAbsoluteMassDeviation
             settingsObject[["FragmentPeakMatchRelativeMassDeviation"]] <- FragmentPeakMatchRelativeMassDeviation
             settingsObject[["MetFragDatabaseType"]] <- MetFragDatabaseType

             settingsObject[["NeutralPrecursorMass"]] <- NeutralPrecursorMass

             settingsObject[["PeakList"]]<- PeakList

             # run MetFrag
             scored.candidates <- metfRag::run.metfrag(settingsObject)

           }
)


#' @title GenerateMetFragLocalDB
#' @author Zhiwei Zhou
#' @description generate candidate list for MetFrag
#' @param metfrag_candidate the result generated from MetFrag
#' @export
#' @examples

setGeneric(name = 'GenerateMetFragLocalDB',
           def = function(
             metfrag_candidate
           ){
             inchikeys <- ImmsTools::SplitInchiKey(inchikey = metfrag_candidate$InChIKey)

             local_db <- data.frame(Identifier=metfrag_candidate$Identifier,
                                    InChI=metfrag_candidate$InChI,
                                    MonoisotopicMass=metfrag_candidate$MonoisotopicMass,
                                    MolecularFormula=metfrag_candidate$MolecularFormula,
                                    InChIKey1=inchikeys$inchikey1,
                                    InChIKey2=inchikeys$inchikey2,
                                    InChIKey3=inchikeys$inchikey3,
                                    InChI=inchikeys$inchikey,
                                    SMILES=metfrag_candidate$SMILES,
                                    Name=metfrag_candidate$IUPACName,
                                    stringsAsFactors = F)

             return(local_db)
           }
)


#' @title GenerateCfmIdCandidates
#' @author Zhiwei Zhou
#' @description generate candidate list for CFM-ID
#' @param metfrag_candidate the result generated from MetFrag
#' @export
#' @examples

setGeneric(name = 'GenerateCfmIdCandidates',
           def = function(
             metfrag_candidate
           ){
             inchikeys <- ImmsTools::SplitInchiKey(inchikey = metfrag_candidate$InChIKey)

             local_db <- data.frame(Identifier=metfrag_candidate$Identifier,
                                    SMILES=metfrag_candidate$SMILES,
                                    stringsAsFactors = F)

             return(local_db)
           }
)


#' @title GenerateCfmIdSpec
#' @author Zhiwei Zhou
#' @description generate peak list for CFM-ID
#' @param energy0 10v specttra.  1st column: mz; 2nd column: intensity.
#' @param energy1 20v specttra.  1st column: mz; 2nd column: intensity.
#' @param energy2 40v specttra.  1st column: mz; 2nd column: intensity.
#' @param dir_path Default: '.'
#' @param file_name Default:'temp_spec.txt'
#' @export
#' @examples
#' GenerateCfmIdSpec(energy)

setGeneric(name = 'GenerateCfmIdSpec',
           def = function(
             energy0=NULL,
             energy1=NULL,
             energy2=NULL,
             dir_path='.',
             file_name='temp_spec.txt'
           ){

             fileConn <- file(description = file.path(dir_path, file_name), open = "a")

             if (!is.null(energy0)) {
               cat(paste('energy0', '\n', sep=''), file = fileConn)
               for (i in 1:nrow(energy0)) {
                 cat(paste(energy0[i, 1], ' ', energy0[i, 2], '\n', sep=''), file = fileConn)
               }
             }

             if (!is.null(energy1)) {
               cat(paste('energy1', '\n', sep=''), file = fileConn)
               for (i in 1:nrow(energy1)) {
                 cat(paste(energy1[i, 1], ' ', energy1[i, 2], '\n', sep=''), file = fileConn)
               }
             }

             if (!is.null(energy2)) {
               cat(paste('energy2', '\n', sep=''), file = fileConn)
               for (i in 1:nrow(energy2)) {
                 cat(paste(energy2[i, 1], ' ', energy2[i, 2], '\n', sep=''), file = fileConn)
               }
             }

             close(fileConn)
           }
)


# cfm-id.exe <spectrum_file> <id> <candidate_file> <num_highest> <ppm_mass_tol> <abs_mass_tol>
#   <prob_thresh> <param_file> <config_file> <score_type> <apply_postprocessing> <output_file> <output_msp_or_mgf>

# cfm-id.exe example_spec.txt TEST_01 example_candidates.txt -1 10.0 0.01 0.001 metab_se_cfm/param_output0.log metab_se_cfm/param_config.txt Jaccard 1 TEST_01.txt TEST_01.msp

#' @title RunCFMID
#' @author Zhiwei Zhou
#' @description R interfaces to run CFM-ID.exe. Detail: https://sourceforge.net/p/cfm-id/wiki/Home/#cfm-id-precomputed
#' @param cfm_id the path of cfm_id.exe Default: 'F:/MetIMMS_MetFrag/cfm_id/cfm-id.exe'
#' @param spec_file the path of spec_file txt file. 1st column: mz; 2nd column: intensity. It could be generated by GenerateCfmIdSpec
#' @param spec_file the path of spec_file txt file. 1st column: identifier; 2nd column: SMILES. It could be generated by GenerateCfmIdSpec
#' @param id An identifier for the target molecule (Used to retrieve input spectrum from msp (if used). Default: 'temp'
#' @param num_highest The number of (ranked) candidates to return or -1 for all. Default: -1
#' @param ppm_mass_tol The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10ppm). Default: 10
#' @param abs_mass_tol The mass tolerance in abs Da to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs ( if not given defaults to 0.01Da). Default: 0.01
#' @param prob_thresh The probability below which to prune unlikely fragmentations (default 0.001). Default: 0.001
#' @param param_file This file is the output of cfm-train. Pre-trained models as used in the above publication can be found in the supplementary data for that paper stored within the source tree of this project. Here, we use the se sfm model for both positive and negative modes. Positive mode: F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_output0.log; Negative mode: F:/MetIMMS_MetFrag/cfm_id/negative_metab_se_cfm/param_output0.log. Default:  'F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_output0.log'
#' @param config_file The filename where the configuration parameters of the cfm model can be found. This needs to match the file passed to cfm-train during training. See cfm-train documentation below for further details. Here, we use the se sfm model for both positive and negative modes. Positive mode: F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_config.txt; Negative mode: F:/MetIMMS_MetFrag/cfm_id/negative_metab_se_cfm/param_config.txt. Default:  'F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_config.txt'
#' @param score_type The type of scoring function to use when comparing spectra. Options: Jaccard (default), DotProduct. Default: Jaccard
#' @param apply_postprocessing Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first) (0 = OFF (default for EI-MS), 1 = ON (default for ESI-MS/MS)). Default: 1
#' @param output_file the path of output result. It should be in ".txt" format. Default: 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_result.txt'
#' @param output_msp the path of output MSP file (MS/MS file). It should be in ".msp" format. Default: 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_spec.msp'
#' @export
#' @examples
#' RunCFMID(spec_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_test_result/example_spec.txt',
#'          candidate_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_test_result/example_candidates.txt',
#'          ppm_mass_tol = 10,
#'          abs_mass_tol = 0.01,
#'          score_type = 'Jaccard',
#'          output_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_result.txt',
#'          output_msp = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_spec.msp')

setGeneric(name = 'RunCFMID',
           def = function(
             cfm_id = 'F:/MetIMMS_MetFrag/cfm_id/cfm-id.exe',
             spec_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_test_result/example_spec.txt',
             id = 'temp',
             candidate_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_test_result/example_candidates.txt',
             num_highest = -1,
             ppm_mass_tol = 10,
             abs_mass_tol = 0.01,
             prob_thresh = 0.001,
             param_file = 'F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_output0.log',
             config_file = 'F:/MetIMMS_MetFrag/cfm_id/metab_se_cfm/param_config.txt',
             score_type = c('Jaccard', 'DotProduct'),
             apply_postprocessing = 1,
             output_file = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_result.txt',
             output_msp = 'F:/MetIMMS_MetFrag/20190521_retrieve_candidates/test/cfm_id_output/TEST_01_spec.msp'
           ){
             score_type <- match.arg(score_type)

             if (!stringr::str_detect(param_file, 'param_output0.log')) {stop('Please check param_file [param_output0.log]\n')}
             if (!stringr::str_detect(config_file, 'param_config.txt')) {stop('Please check config_file [param_config.txt]\n')}
             if (!stringr::str_detect(output_msp, '.msp')) {stop('Please check the file name of output msp\n')}


             cmdl <- paste(cfm_id, spec_file, id, candidate_file, num_highest, ppm_mass_tol,
                           abs_mass_tol, prob_thresh, param_file, config_file, score_type, apply_postprocessing,
                           output_file, output_msp, sep=' ')

             shell(cmdl, translate = TRUE)

             cat('CFM-ID has completed\n')
           }
)
