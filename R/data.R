#' Atmar & Patterson (1995) incidence matrices
#'
#' This is a list of site by species binary *incidence* (0 or 1) matrices.
#' Note that after all cleaning procedures, the final number of matrices is 290 (not 294).
#' Each matrix in the list has a name that enables it to be linked back to the
#' original dataset in the publication.
#'
#' @source Atmar W. & Patterson B. D. (1995) The nestedness temperature calculator:
#' a Visual Basic  program,  including  294  presence absence matrices. AICS Research,
#' Univ. Park, NM and Field Museum, Chicago, http://aics-research.com/nestedness/.
"data.Atmar"


#' Ulrich & Gotelli (2010) abundance matrices
#'
#' This is a list of 186 site by species *abundance* matrices.
#' Each matrix in the list has a name that enables it to be linked back to the
#' original dataset in the publication.
#'
#' @source Ulrich W. & Gotelli N.J. (2010) Null model analysis of species associations
#' using abundance data. Ecology 91: 3384-3397.
"data.Ulrich"


#' Orwig, Foster & Ellison (2015) Harvard forest plot
#'
#' This is a set of all 'alive' trees from the file hf253-03-trees-2014.csv, downloaded from
#' http://harvardforest.fas.harvard.edu:8080/exist/apps/datasets/showData.html?id=hf253
#' This represents the 2014 census.
#'
#' @source Owig D., Foster D. & Ellison A. (2015) Harvard Forest CTFS-ForestGEO
#' Mapped Forest Plot since 2014. Harvard Forest Data Archive: HF253.
"data.Harvard"

#' Number of papers on biodiversity + a given category on Clarivate Web of Science
#'
#' The search was done on 19 Sep 2019. The search categories were the given term (in the data)
#' AND "biodiversity".
#'
#' @source http://apps.webofknowledge.com/
"data.WOS"

#' Data on 3879 papers published in American Nautralist, Ecology, and Ecography
#' in years 1995, 1999, 2003, 2007, 2011, and 2019
#'
#' I created the data by
#' reading 3879 abstracts and manually classifying them into thematic categories.
#' The abstracts were taken from Clarivate Web of Science. Metadata are in the "data" folder.
#'
#' @source http://apps.webofknowledge.com/
"data.journals"
