allDuplicatedEntries <- function(vektor) 
{
  if (length(vektor) == 0) 
    return(0)
  vektab = data.table::data.table(myvektor = vektor, num = 1:length(vektor))
  duplicated_vals = vektab[duplicated(myvektor), myvektor]
  duplicated_entries = vektab[myvektor %in% duplicated_vals]
  data.table::setkey(duplicated_entries, myvektor)
  duplicated_entries$num
}
