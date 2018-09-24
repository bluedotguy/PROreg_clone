\name{COPD}
\alias{COPD}
\title{COPD: chronic obstructive pulmonary disease}
\usage{
data(COPD)
}
\format{The \code{COPD} data frame contains the recoded scores based on
\describe{
\item{\code{RP.}}{Role Physical. Integer values [0,1,2,3,4]}

\item{\code{GH.}}{General Health. Integer values [0,1,2,...,20]}

\item{\code{RE.}}{Role Emotional. Integer values [0,1,2,3]}

\item{\code{time}}{Time points. Integer values [1,2,3,4]}

\item{\code{sex}}{sex 1=male; 2=female}

\item{\code{id}}{Individual, total of 543}

\item{cluster}{cluster A,B,C,D}

\item{dd}{Days (integer)}

\item{dy}{Days/365}
}
}
\details{
COPD data from Najera-Zuloaga et al. (2018).
}
\description{
The COPD Study is an observational study that was designed at the Respiratory Unit at Galdakao Hospital in Spain. Patients' COPD had to be stable for at least six months before enrolment, and they were followed up to five years with clinical examination and interviews. The health-status was measured with the generic questionnaire SF-36 Health Survey. The \code{COPD} data frame contains three dimensions (Role Physical, General Health and Role Emotional) of 543 patients recoded as a binomial form according to Arostegui et al. (2013). See function \code{?SF36rec}. Patients are classified into 4 clusters or subtypes: A,B, C and D. Subtypes A, B, and C, had marked respiratory profiles with a continuum in severity of several variables and subtype D, had a more systemic profile with intermediate respiratory disease severity (see Esteban et al. 2016).
}
\examples{
data(COPD)
head(COPD)
str(COPD)
}
\references{
\itemize{
   \item Arostegui, I., Nunez-Anton, V., and Quintana, J. M. (2013). "On the recoding of continuous and bounded indexes to a binomial form: an application to quality-of-life scores". Journal of Applied Statistics, 40:563-582. \href{https://doi.org/10.1080/02664763.2012.749845}{https://doi.org/10.1080/02664763.2012.749845}
   \item Esteban, C., Arostegui, I., Aburto, M., Moraza, J., Quintana, J. M., Garcia-Loizaga, A., Basualdo, L. V., Aramburu, A., Aizpiri, S., Uranga, A., and Capelastegui, A. (2016). \href{https://doi.org/10.1371/journal.pone.0161710}{"Chronic obstructive pulmonary disease subtypes. Transitions over time"}. PLoS One, 11(9):e0161710.
   \item Najera-Zuloaga, J., Lee, D.-J., and Arostegui, I. (2018). \href{https://doi.org/10.1177/0962280217690413}{"Comparison of beta-binomial regression model approaches to analyze health-related quality of life data"}. Stat Methods Med Res. 2018 Oct:27(10):2989-3009.
   \item Najera-Zuloaga J., Lee D.-J. & Arostegui I. (2018): A beta-binomial mixed-effects model approach for analysing longitudinal discrete and bounded outcomes, to appear in \emph{Biometrical Journal}.
 }
}
