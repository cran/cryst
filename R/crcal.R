#' @title Relative Crystallinity Calculation of X-Ray Diffraction Pattern of Starch by Bruckner Method
#' @details Calculate the relative starch crystallinity of XRD pattern by Bruckner method.
#' @description Allow  to calculate the relative crystallinity of starch by XRD.
#' The basic concept of Bruckner approach involves obtaining a smoothed
#' line that separates the amorphous and crystalline starch of an X-ray
#' diffraction pattern. This smoothed line is achieved by applying
#' a moving average smoothing method to the original pattern.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param pattern matrix. The matrix of X-ray diffraction pattern.  The first row corresponds to Bragg angle 2\eqn{\theta}; the second row corresponds to intensity.
#' @param N numeric. N length of the smoothing window (number of variables). Defaults to 41.
#' @param iter numeric. Iter number of iterations. Defaults to 100.
#' @return An object of class crystMW, which is a list with the following components:
#'   \item{original}{Original matrix of X-ray diffraction patterns.}
#'   \item{background}{Estimation of the background shape (curve of the amorphous starch).}
#'   \item{corrected}{Estimation of residual crystalline area (curve of the crystalline starch).}
#'   \item{summary}{Summary calculation of crystallinity. Total area under the curve of the diffraction pattern (A.U.); Amorphous area (A.U.); Crystalline area (A.U.); Relative crystallinity (\%).}
#' @importFrom flux auc
#' @examples
#' data(XRD)
#' # Convert data frame to matrix, select A-type starch
#' pattern <- as.matrix(t(XRD[,c("Bragg_angle","A")]))
#' # List of crystallinity components
#' crs <- crystMW(pattern, N = 41, iter = 100)
#' # Original matrix
#' original <- crs$original
#' # Background shape
#' background <- crs$background
#' # Curve of the crystalline starch
#' corrected <- crs$corrected
#' # Summary calculation of crystallinity
#' summary <- crs$summary
#' @references
#' Bruckner, S. (2000). "Estimation of the background in powder diffraction patterns through a robust smoothing procedure." Journal of Applied Crystallography 33(3 Part 2): 977-979.
#' @export
crystMW<-function(pattern, N=41, iter=100){
  if (pattern[1, 1] < pattern[1, dim(pattern)[2]]) {
    pattern <- t(apply(pattern, 1, rev))
  }
  if (missing(pattern)) {
    stop('No pattern data provided')
  }
  if (missing(N)) {
    stop('No length of the smoothing window provided')
  }
  if (missing(iter)) {
    stop('No number of iterations provided')
  }
  # Inserting on both sides N points
  s  <- pattern[2,]
  sa <- rep(NA, length(s)+2*N)
  sa <- c(rep(pattern[2,1],N),s,rep(pattern[2,length(s)],N))
  # Iteration
  for(it in 1:iter){
    # Replacement of the kth point by the average of its 2N next neighbours (N to the left, N to the right)
    Ik <- sa
    for(k in N:(length(sa)-N)){
      Ik[k] <- (1/(2*N))*(sum(sa[(k-N):(k-1)])+sum(sa[(k+1):(k+N)]))
    }
    # The new pattern is compared, point by point, with the original one, retaining the lowest of the two intensities
    for(i in 1:length(sa)){
      if(Ik[i] <= sa[i]){
        sa[i] <- Ik[i]
      }else{
        sa[i] <- sa[i]
      }
    }
  }
  # Reconstructed background
  background <-as.matrix(rbind(pattern[1,], sa[(N+1):(length(sa)-N)]))
  colnames(background) <- NULL
  rownames(background) <- NULL
  # Corrected pattern (original minus background)
  corrected<-as.matrix(rbind(pattern[1,], s-background[2,]))
  colnames(corrected) <- NULL
  rownames(corrected) <- NULL
  # Crystalline area
  res <- matrix(NA, nrow=1, ncol=4)
  colnames(res) <-c('Total Area', 'Amorphous area', 'Crystalline area', 'Relative crystallinity')
  rownames(res) <- NULL
  ## Total area
  res[1] <- auc(pattern[1,], pattern[2,], thresh = 0, dens = 100)
  ## Amorphous area
  res[2] <- auc(background[1,], background[2,], thresh = 0, dens = 100)
  ## Crystalline area
  res[3] <-res[1]-res[2]
  ## Relative crystallinity
  res[4] <- res[3]/res[1]
  # Output
  value <-list(original=pattern, background=background, corrected=corrected, summary=res)
  attr(value, "class") <- "crystMW"
  value
}

#' @title Relative Crystallinity Calculation of X-Ray Diffraction Pattern of Starch by Frost Method
#' @details Calculate the relative starch crystallinity of XRD pattern by Frost method.
#' @description Allow  to calculate the relative crystallinity of starch by XRD.
#' The basic concept of Frost approach involves obtaining a smoothed
#' line that separates the amorphous and crystalline starch of an X-ray
#' diffraction pattern. This smoothed line is achieved by applying
#' a Savitzky-Golay smoothing method to the original pattern.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param pattern matrix. The matrix of X-ray diffraction pattern.  The first row corresponds to Bragg angle 2\eqn{\theta}; the second row corresponds to intensity.
#' @param N numeric. N length of the smoothing window (number of variables). Defaults to 101.
#' @param iter numeric. Iter number of iterations. Defaults to 400.
#' @param p numeric. Filter order. Defaults to 2.
#' @return An object of class crystSG, which is a list with the following components:
#'   \item{original}{Original matrix of X-ray diffraction patterns.}
#'   \item{background}{Estimation of the background shape (curve of the amorphous starch).}
#'   \item{corrected}{Estimation of residual crystalline area (curve of the crystalline starch).}
#'   \item{summary}{Summary calculation of crystallinity. Total area under the curve of the diffraction pattern (A.U.); Amorphous area (A.U.); Crystalline area (A.U.); Relative crystallinity (\%).}
#' @importFrom flux auc
#' @importFrom stats convolve
#' @importFrom pracma pinv
#' @examples
#' data(XRD)
#' # Convert data frame to matrix, select A-type starch
#' pattern <- as.matrix(t(XRD[, c("Bragg_angle","A")]))
#' # List of crystallinity components
#' crs <- crystSG(pattern, N = 101, iter = 400, p = 2)
#' # Original matrix
#' original <- crs$original
#' # Background shape
#' background <- crs$background
#' # Curve of the crystalline starch
#' corrected <- crs$corrected
#' # Summary calculation of crystallinity
#' summary <- crs$summary
#' @references
#' Frost, K., et al. (2009). "Crystallinity and structure of starch using wide angle X-ray scattering." Carbohydrate Polymers 78(3): 543-548.
#' @export
crystSG<-function(pattern, N=101, iter=400, p=2){
  if (pattern[1, 1] < pattern[1, dim(pattern)[2]]) {
    pattern <- t(apply(pattern, 1, rev))
  }
  if (missing(pattern)) {
    stop('No pattern data provided')
  }
  if (missing(N)) {
    stop('No length of the smoothing window provided')
  }
  if (missing(iter)) {
    stop('No number of iterations provided')
  }
  if (missing(p)) {
    stop('No filter order provided')
  }
  # Inserting on both sides N points
  s  <- pattern[2,]
  sa <- rep(NA, length(s)+2*N)
  sa <- c(rep(pattern[2,1],N),s,rep(pattern[2,length(s)],N))
  # Iteration
  for(it in 1:iter){
    # Replacement of the kth point by the savitzky golay smoothing of its 2N next neighbours (N to the left, N to the right)
    fc <- (N - 1)/2
    X  <- outer(-fc:fc, 0:p, FUN = "^")
    Y  <- pinv(X)
    t  <- convolve(sa, rev(Y[1,]), type = "o")
    Ik <- t[(fc + 1):(length(t) - fc)]
    # The new pattern is compared, point by point, with the original one, retaining the lowest of the two intensities
    for(i in 1:length(sa)){
      if(Ik[i]<=sa[i]){
        sa[i]<-Ik[i]
      }else{
        sa[i]<-sa[i]
      }
    }
  }
  # Reconstructed background
  background <-as.matrix(rbind(pattern[1,], sa[(N+1):(length(sa)-N)]))
  colnames(background) <- NULL
  rownames(pattern)<- NULL
  # Corrected pattern (original minus background)
  corrected<-as.matrix(rbind(pattern[1,], s-background[2,]))
  colnames(corrected) <- NULL
  rownames(corrected) <- NULL
  # Crystalline area
  res <- matrix(NA, nrow=1, ncol=4)
  colnames(res) <-c('Total Area', 'Amorphous area', 'Crystalline area', 'Relative crystallinity')
  rownames(res) <- NULL
  ## Total area
  res[1] <- auc(pattern[1,], pattern[2,], thresh = 0, dens = 100)
  ## Amorphous area
  res[2] <- auc(background[1,], background[2,], thresh = 0, dens = 100)
  ## Crystalline area
  res[3] <-res[1]-res[2]
  ## Relative crystallinity
  res[4] <- res[3]/res[1]
  # Output
  value <-list(original=pattern, background=background, corrected=corrected, summary=res)
  attr(value, "class") <- "crystSG"
  value
}

#' @title Relative Crystallinity Calculation of FTIR Spectrum of Starch by SUN Method
#' @details Calculate the relative starch crystallinity of FTIR spectrum by SUN method.
#' @description Allow  to calculate the relative crystallinity of starch by FTIR.
#' The basic concept of SUN approach involves obtaining a gaussian holocrystalline-peak
#' in the 800-1300 cm-1 region of FTIR spectrum of starch which is divided into amorphous
#' region and crystalline region.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectrum matrix. The matrix of FTIR spectrum baseline-corrected by drawing a tangentline in the 800-1300 cm-1 region. The first row corresponds to wavelength; the second row corresponds to intensity.
#' @param mu numeric. Gaussian mean of holocrystalline-peak. Defaults to 1180.
#' @param sigma numeric. Standard deviation of holocrystalline-peak. Defaults to 60.
#' @param k numeric. Arbitrary scaling parameter. Defaults to 1.
#' @param lim vector. Fitting points of holocrystalline-peak. Defaults to c(1190, 1160, 985, 950).
#' @return An object of class fitFTIRc, which is a list with the following components:
#'   \item{original}{Original matrix of FTIR spectrum.}
#'   \item{gauss}{Gaussian curve fit.}
#'   \item{fit}{Summary of Non-Linear Least-Squares Model Fits.}
#'   \item{summary}{Summary calculation of crystallinity. Total area under the curve of the diffraction spectrum (A.U.); Amorphous area (A.U.); Crystalline area (A.U.); Relative crystallinity (\%).}
#' @importFrom flux auc
#' @importFrom stats nls
#' @importFrom stats coef
#' @examples
#' # Convert data frame to matrix, select A-type starch
#' spectrum <- as.matrix(t(FTIR[, c('wavelength','A')]))
#' # List of crystallinity components
#' crs <- fitFTIRc(spectrum = spectrum, mu = 1180, sigma = 60, k = 1, lim = c(1190, 1160, 985, 955))
#' # Original matrix
#' original <- crs$original
#' # Gaussian curve fit
#' gauss <- crs$gauss
#' # Summary of Non-Linear Least-Squares Model Fits
#' fit <- crs$fit
#' # Summary calculation of crystallinity
#' summary <- crs$summary
#' @references
#' Sun, Y., et al. (2014). "A new method for determining the relative crystallinity of chickpea starch by Fourier-transform infrared spectroscopy." Carbohydrate Polymers 108: 153-158.
#' @export
fitFTIRc<- function(spectrum, mu=1180, sigma=60, k=1, lim=c(1190, 1160, 985, 950)){
  if (spectrum[1, 1] < spectrum[1, dim(spectrum)[2]]) {
    spectrum <- t(apply(spectrum, 1, rev))
  }
  if (missing(spectrum)) {
    stop('No spectrum data provided')
  }
  if (missing(mu)) {
    stop('No gaussian mean provided')
  }
  if (missing(sigma)) {
    stop('No standard deviation provided')
  }
  if (missing(k)) {
    stop('No arbitrary scaling parameter provided')
  }
  if (missing(lim)) {
    stop('No fitting points of holocrystalline-peak provided')
  }
  # Truncated region between 950-985 cm-1 and 1160-1190 cm-1
  lim<-sort(lim, decreasing = TRUE)
  li01 <- which(round(spectrum[1, ], 0) == lim[1])[1]
  ls01 <- which(round(spectrum[1, ], 0) == lim[2])[1]
  li02 <- which(round(spectrum[1, ], 0) == lim[3])[1]
  ls02 <- which(round(spectrum[1, ], 0) == lim[4])[1]
  a01 <- spectrum[, ls01:li01]
  a02 <- spectrum[, ls02:li02]
  a03 <- cbind(c(1300,rep(0, time=1)),a02,a01, c(800,rep(0, time=1)))
  #Apply function nls
  ## First present the data in a data frame
  s <- data.frame(x=a03[1,], r=a03[2,])
  fit <- nls(r ~ (k*exp(-(x-mu)**2/(2*sigma**2))),data=s, start=list(k=2, mu=1063, sigma=50))
  # Crystalline area
  res <- matrix(NA, nrow=1, ncol=4)
  colnames(res) <-c('Total Area', 'Amorphous area', 'Crystalline area', 'Relative crystallinity')
  rownames(res) <- NULL
  ## Gaussian line polygon
  gl <- matrix(NA, nrow=2, ncol=ncol(spectrum))
  gl[1,] <- spectrum[1,]
  for (i in 1:length(spectrum[1,])){
    gl[2,i]<-coef(fit)[1] * exp(-(gl[1,i]-coef(fit)[2])**2/(2*coef(fit)[3]**2))
  }
  ## Crystalline polygon
  pol <- matrix(NA, nrow=2, ncol=ncol(spectrum))
  pol[1,] <- spectrum[1,]
  for (i in 1:length(spectrum[1,])){
    if (gl[2,i]<=spectrum[2,i]){
      pol[2,i]<-gl[2,i]
    }else{
      pol[2,i]<-spectrum[2,i]
    }
  }
  ## Total area
  res[1] <- auc(gl[1,], gl[2,], thresh = 0, dens = 100) # Gaussian line polygon
  ## Crystalline area
  res[2] <- auc(pol[1,], pol[2,], thresh = 0, dens = 100) # Crystalline polygon
  ## Amorphous area
  res[3] <-res[1]-res[2]
  ## Relative crystallinity
  res[4] <- res[2]/res[1]
  # Output
  value <- list(original=spectrum, gauss=gl, fit=summary(fit), summary=res)
  attr(value, "class") <- "fitFTIRc"
  value
}

#' @title Plots the Crystalline Area of a X-Ray Diffraction Pattern of Starch
#' @description Produces a graph of the crystalline area of a X-ray diffraction pattern of starch and background.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param pattern matrix. The matrix of X-ray diffraction pattern.  The first row corresponds to Bragg angle 2\eqn{\theta}; the second row corresponds to intensity.
#' @param background matrix. The matrix of background shape (curve of the amorphous starch). The first row corresponds to Bragg angle 2\eqn{\theta}; the second row corresponds to intensity.
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics polygon
#' @importFrom graphics title
#' @importFrom grDevices rgb
#' @examples
#' # Convert data frame to matrix, select A-type starch
#' pattern <- as.matrix(t(XRD[, c("Bragg_angle","A")]))
#' # List of crystallinity components
#' crs <- crystMW(pattern, N = 11, iter = 100)
#' # Original matrix
#' original <- crs$original
#' # Background shape
#' background <- crs$background
#' # Plots the crystalline area of a XRD pattern
#' xrdplot(pattern=original, background=background)
#' @export
xrdplot <- function(pattern, background){
  if (missing(pattern)) {
    stop('No pattern data provided')
  }
  if (missing(background)) {
    stop('No background provided')
  }
  if (pattern[1, 1] < pattern[1, dim(pattern)[2]]) {
    pattern <- t(apply(pattern, 1, rev))
  }
  if (background[1, 1] < background[1, dim(background)[2]]) {
    background <- t(apply(background, 1, rev))
  }

  titx= '2'~{theta}         # x title
  tity= 'Intensity (a.u.)'  # y title
  cexLab= 1.5

  # Limits
  xlimS <- round(max(pattern[1,]), 0)
  xlimI <- round(min(pattern[1,]), 0)
  ylimS <- round(max(pattern[2, ])*1.3,1)
  ylimI <- round(min(background[2, ]),2)
  xlim<-c(xlimI, xlimS)
  ylim<-c(ylimI, ylimS)

  # Plot XRD
  plotpattern<-plot(pattern[1, ], pattern[2, ],  # function to plot
                    xlab="",          # suppress x labels
                    ylab="",          # suppress y labels
                    type = 'l',       # specify line graph
                    xlim = xlim,      # extend axis limits to give space for text annotation
                    ylim = ylim       # ditto
  )
  lines(background[1, ], background[2, ],col="black")
  polygon(c(pattern[1,],rev(background[1,])),c(pattern[2,],rev(background[2,])),col="grey")

  # Titles
  title(xlab= titx, col.lab=rgb(0, 0, 0),  cex.lab=cexLab, line = 3.8)
  title(ylab= tity, col.lab=rgb(0, 0, 0),  cex.lab=cexLab)
}

#' @title Plots the Crystalline Area of a FTIR Spectrum of Starch
#' @description Produces a graph of the crystalline area of a FTIR spectrum of starch and the Gauss curve.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectrum matrix. The matrix of FTIR spectrum baseline-corrected by drawing a tangentline in the 800-1300 cm-1 region. The first row corresponds to wavelength; the second row corresponds to intensity.
#' @param gauss matrix. The matrix of Gauss curve (gaussian holocrystalline-peak).
#' @param lim vector. Regions of the FTIR spectrum comprising the fixing points of the Gauss curve.
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics abline
#' @importFrom graphics title
#' @importFrom graphics polygon
#' @importFrom grDevices rgb
#' @examples
#' # Convert data frame to matrix, select A-type starch
#' spectrum <- as.matrix(t(FTIR[, c('wavelength','A')]))
#' # List of crystallinity components
#' crs <- fitFTIRc(spectrum = spectrum, mu = 1180, sigma = 60, k = 1, lim = c(1190, 1160, 985, 955))
#' # Original matrix
#' original <- crs$original
#' # Gaussian curve fit
#' gauss <- crs$gauss
#' # Plots the crystalline area of a FTIR spectrum
#' ftirplot(spectrum=original, gauss=gauss, lim=c(1190, 1160, 985, 950))
#' @export
ftirplot <- function(spectrum, gauss, lim=c(1190, 1160, 985, 950)) {
  if (missing(spectrum)) {
    stop('No spectrum data provided')
  }
  if (missing(gauss)) {
    stop('No gauss curve provided')
  }
  if (missing(lim)) {
    stop('No fitting points of holocrystalline-peak provided')
  }
  if (spectrum[1, 1] < spectrum[1, dim(spectrum)[2]]) {
    spectrum <- t(apply(spectrum, 1, rev))
  }
  if (gauss[1, 1] < gauss[1, dim(gauss)[2]]) {
    gauss <- t(apply(gauss, 1, rev))
  }
  # Titles
  titx= 'Wavelength'~(cm^{-1}) # x title
  tity= "Absorbance (a.u.)"    # y title
  cexLab= 1.5                  # magnification of x and y labels
  # Limits
  xlimS <- round(max(spectrum[1,]), -1)
  xlimI <- round(min(spectrum[1,]), -1)
  ylimS <- round(max(gauss[2,])*1.2,1)
  ylimI <- round(min(spectrum[2, ]),2)
  xlim<-c(xlimS, xlimI)
  ylim<-c(ylimI, ylimS)
  # Plot spectrum and abline
  plot(gauss[1,], gauss[2,], type='l',xlab="", ylab="",yaxt='n', xlim=xlim, ylim=ylim)
  lines(spectrum[1,],spectrum[2,],col="black")
  lim<-sort(lim, decreasing = TRUE)
  abline(v=lim[4], lty= 2)
  abline(v=lim[3], lty= 2)
  abline(v=lim[2], lty= 2)
  abline(v=lim[1], lty= 2)
  # Calculate polygon of crystalline starch
  pol <- matrix(NA, nrow=2, ncol=ncol(spectrum))
  pol[1,] <- spectrum[1,]
  for (i in 1:length(spectrum[1,])){
    if (gauss[2,i]<=spectrum[2,i]){
      pol[2,i]<-gauss[2,i]
    }else{
      pol[2,i]<-spectrum[2,i]
    }
  }
  # Crystalline polygon
  pol <- matrix(NA, nrow=2, ncol=ncol(spectrum))
  pol[1,] <- spectrum[1,]
  for (i in 1:length(spectrum[1,])){
    if (gauss[2,i]<=spectrum[2,i]){
      pol[2,i]<-gauss[2,i]
    }else{
      pol[2,i]<-spectrum[2,i]
    }
  }
  # Titles
  title(xlab= titx, col.lab=rgb(0, 0, 0),  cex.lab=cexLab, line = 3.4)
  title(ylab= tity, col.lab=rgb(0, 0, 0),  cex.lab=cexLab, line = 1)
  # Apply cross-hatching to a polygon
  polygon(pol[1, ], pol[2, ], density=10, angle=45, border=NULL)
}
