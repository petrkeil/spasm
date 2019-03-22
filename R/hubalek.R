#' Hubalek's indices calculated all at once
#'
#' This is a function that calculates all of the indices listed in Hubalek (1982).
#'
#' @param x1 A binary vector.
#' @param x2 A binary vector of the same legnth as x1.
#'
#' @return A dist or data.frame objects with the pairwise association values. Letters and
#' numbers correspond to those in Hubalek's paper.
#' @references Hubalek Z. (1982) Coefficients of association and similarity,
#' based on binary (presence-absence) data: and evaluation. Biol. Rev., 57: 669-689.
#' @export
#'

hubalek <- function(x1, x2)
{

  A = sum((x1+x2 == 2))
  B = sum((x1+x2 == 1)*x1)
  C = sum((x1+x2 == 1)*x2)
  D = sum((x1+x2 == 0))
  N = length(x1)

  res <- numeric(43)
  names(res) <- paste("A", 1:43, sep="")

  # the actual indices:
  res[1]  = A/(max(A + B, A + C)) # Braun-Blanquet (1932)
  res[2]  = A/(min(A + B, A + C)) # Simpson (1943)
  res[3]  = A/(B + C) # Kulczynski (1927)
  res[4]  = A/(A + B + C) # Jaccard (1901)
  res[5]  = A/(A + 0.5*(B+C)) # Dice (1945), Sorensen (1948)
  res[6]  = A/(A + 2*(B + C)) # Sokal & Sneath (1963)
  res[7]  = 0.5*(A/(A+B) + A/(A+C)) # Kulczynski (1927)
  res[8]  = (A/2)*(1/(A+B) + 1/(A+C)) # Driver & Kroeber (1932)
  res[9]  = A/(A+B) + A/(A+C)# Johnson (1967)
  res[10] = (A*A - B*C)/((A+B)*(A+C)) # McConnaughey (1964)
  res[11] = A / (((A+B)*(A+C))^0.5 ) # Driver & Kroeber (1932), Ochiai (1957)
  res[12] = (A*A)/((A+B)*(A+C)) # Sorgenfrei (1959)
  res[13] = A/(((A+B)*(A+C))^0.5) - 0.5*(max(A+B, A+C)) # Fager & McGowan (1963)
  res[14] = A/(A + B + C + D) # Russell & Rao (1940)
  res[15] = A/(0.5*(A*B + A*C) + B*C) # Mountford (1962)
  res[16] = (A*B + B*C) / (A*B + 2*B*C + C*D) # Peirce (1884)
  res[17] = (A - (A+B)*(A+C)) / ((A+B)*(C+D)*(A+C)*(B+D)) # Eyraud (1936)
  res[18] = 0.25*(A/(A+B) + A/(A+C) + D/(C+D) + D/(B+D)) # Sokal & Sneath (1963)
  res[19] = (A+D) / (B+C) # Sokal & Sneath (1963)
  res[20] = (A + D) / N # Sokal & Michener (1958)
  res[21] = (50*pi)^(-1)* asin(((A + D) / N)^0.5) # Goodall (1967), Austin & Colwell (1977)
  res[22] = (A+D)/(A + 0.5*(B+C) + D) # Sokal & Sneath (1963)
  res[23] = (A+D)/(A + 2*(B+C+D)) # Rogers & Tanimoto (1960)
  res[24] = (A+D-B-C)/N # Hamann (1961)
  res[25] = A*D/(((A+B)*(C+D)*(A+C)*(B+D))^0.5) # Sokal & Sneath (1963)
  res[26] = (A*D - B*C)/((A+C)*(B+D)) # Pierce (1884)
            Chi.sq = N*((A*D - B*C)^2)/((A+B)*(C+D)*(A+C)*(B+D)) # Pearson (1905)
  res[27] = Chi.sq
  res[28] = Chi.sq/((N + Chi.sq)^0.5) # Pearson (1905)
  res[29] = sqrt(2) * (A*D - B*C) / (((A*D - B*C)^2 - (A+B)*(C+D)*(A+C)*(B+D))^0.5) # Cole (1949)
  res[30] = (A*D + B*C) / (((A+B)*(C+D)*(A+C)*(B+D))^0.5) # Yule (1912), Pearson & Heron (1913)
  res[31] = (A*D - B*C)^2 / ((A+B)*(C+D)*(A+C)*(B+D)) # Doolittle (1885), Pearson (1926)
  res[32] = (sqrt(A*D) + A)/(sqrt(A*D) + A + B + C) # Baroni-Urbani & Buser (1976)
  res[33] = (sqrt(A*D) + A - B - C)/(sqrt(A*D) + A + B + C)
  res[34] = NA # I don't 100% understand the description here. Cole (1949)
  res[35] = NA # I don't 100% understand the description here. Cole (1949)
  res[36] = (A*D - B*C) / (A*D + B*C) # Yule (1900)
  res[37] = (sqrt(A*D) - sqrt(B*C)) / (sqrt(A*D) + sqrt(B*C)) # Yule (1900)
  res[38] = cos(180*sqrt(B*C) / (sqrt(A*D) + sqrt(B*C))) # Pearson & Heron (1913)
  res[39] = 4*(A*D - B*C)/( (A+D)^2 + (B+C)^2 ) # Michael (1920)
  res[40] = N*A/((A+B)*(A+C)) # Forbes(1907)
  res[41] = log(A) - log(N) - log((A+B)/N) - log((A+C)/N) # Gilbert & Wells (1966)
  res[42] = (N*A -(A+B)*(A+C)) / (N * min(A+B, A+C) - (A+B)*(A+C)) # Forbes (1925)
  res[43] = (N*A -(A+B)*(A+C)) / (N*A + (A+B)*(A+C)) # Tarwid (1960)

  return(res)
}
