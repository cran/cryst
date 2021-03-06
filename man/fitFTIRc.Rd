% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crcal.R
\name{fitFTIRc}
\alias{fitFTIRc}
\title{Relative Crystallinity Calculation of FTIR Spectrum of Starch by SUN Method}
\usage{
fitFTIRc(spectrum, mu = 1180, sigma = 60, k = 1, lim = c(1190, 1160,
  985, 950))
}
\arguments{
\item{spectrum}{matrix. The matrix of FTIR spectrum baseline-corrected by drawing a tangentline in the 800-1300 cm-1 region. The first row corresponds to wavelength; the second row corresponds to intensity.}

\item{mu}{numeric. Gaussian mean of holocrystalline-peak. Defaults to 1180.}

\item{sigma}{numeric. Standard deviation of holocrystalline-peak. Defaults to 60.}

\item{k}{numeric. Arbitrary scaling parameter. Defaults to 1.}

\item{lim}{vector. Fitting points of holocrystalline-peak. Defaults to c(1190, 1160, 985, 950).}
}
\value{
An object of class fitFTIRc, which is a list with the following components:
  \item{original}{Original matrix of FTIR spectrum.}
  \item{gauss}{Gaussian curve fit.}
  \item{fit}{Summary of Non-Linear Least-Squares Model Fits.}
  \item{summary}{Summary calculation of crystallinity. Total area under the curve of the diffraction spectrum (A.U.); Amorphous area (A.U.); Crystalline area (A.U.); Relative crystallinity (\%).}
}
\description{
Allow  to calculate the relative crystallinity of starch by FTIR.
The basic concept of SUN approach involves obtaining a gaussian holocrystalline-peak
in the 800-1300 cm-1 region of FTIR spectrum of starch which is divided into amorphous
region and crystalline region.
}
\details{
Calculate the relative starch crystallinity of FTIR spectrum by SUN method.
}
\examples{
# Convert data frame to matrix, select A-type starch
spectrum <- as.matrix(t(FTIR[, c('wavelength','A')]))
# List of crystallinity components
crs <- fitFTIRc(spectrum = spectrum, mu = 1180, sigma = 60, k = 1, lim = c(1190, 1160, 985, 955))
# Original matrix
original <- crs$original
# Gaussian curve fit
gauss <- crs$gauss
# Summary of Non-Linear Least-Squares Model Fits
fit <- crs$fit
# Summary calculation of crystallinity
summary <- crs$summary
}
\author{
Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
}
\references{
Sun, Y., et al. (2014). "A new method for determining the relative crystallinity of chickpea starch by Fourier-transform infrared spectroscopy." Carbohydrate Polymers 108: 153-158.
}

