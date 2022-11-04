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
  # define variables from params
  if (!is.null(c_start) & keep_fixed) {
    contact_matrix <- params$c_start
    flag_open <- 0
  } else {
    thresh_lockdown <- params$thresh_l
    thresh_open <- params$thresh_o

    c_lockdown <- params$c_lockdown
    c_open <- params$c_open

    # use simpler conditions where measures are only relaxed and not re-tightened
    # for flags
    if (criteria >= thresh_lockdown) {
      flag_open <- 1
    }
    if (criteria <= thresh_open) {
      flag_open <- 2
    }
    
    # for contact matrix
    if (flag_open == 0) {
      contact_matrix <- c_start
    } else if (flag_open == 1) {
      contact_matrix <- c_lockdown
    } else {contact_matric <- c_open}
  }

  rtn <- list(
    contact_matrix = contact_matrix,
    flag_open = flag_open
  )
  return(rtn)
}
