#' Determine contact matrix based on thresholds of cases or IC admissions
#' @param params list of parameter values
#' @param criteria criteria by which to change contact matrix. User can choose between a number of
#' cases per day or a number of IC admissions.
#' @param flag_open an integer value that becomes non-zero (positive) when the contact matrix consistent with
#' "relaxed" non-pharmaceutical interventions is being used.
#' @param keep_fixed logical. if TRUE the contact matrix stays fixed over the
#' entire simulation period
#' @return List of summary results
#' @keywords vacamole
#' @export

choose_contact_matrix_simple <- function(params,criteria,flag_open,keep_fixed) {
  
  # initialise flag and contact matrix
  flag_open <- 0
  contact_matrix <- params$c_start
  
  # define variables from params
  if (keep_fixed) {
    contact_matrix <- params$c_start
  } else if (criteria >= params$thresh_l) {
      flag_open <- 1
      contact_matrix <- params$c_lockdown
  } else if (criteria <= params$thresh_o & flag_open > 0) {
      flag_open <- 2
      contact_matrix <- params$c_open
  }
    
    # for contact matrix
    # if (flag_open == 0) {
    #   contact_matrix <- c_start
    # } else if (flag_open == 1) {
    #   contact_matrix <- params$c_lockdown
    # } else {contact_matric <- params$c_open}
  
  rtn <- list(
    contact_matrix = contact_matrix,
    flag_open = flag_open
  )
  return(rtn)
}
