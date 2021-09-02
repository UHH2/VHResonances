import ROOT
from array import array
# fNorms = {}

# def lnN(nominal, varUp, varDown):
#     if varUp==0: raise Exception("varUp==0")
#     if varDown==0: raise Exception("varDown==0")
#     return (1+abs(nominal-varUp)/nominal,1+abs(nominal-varDown)/nominal)

def lnN(nominal, variation):
    if variation==0: raise Exception("var==0")
    # return 1+abs(nominal-variation)/nominal
    return 1+(nominal-variation)/nominal

def CountBinsMinMax(hist, min, max):
    Nbins = 0
    xmin = 1e6
    xmax = 0
    for i in range(hist.GetNbinsX()+1):
        if hist.GetXaxis().GetBinLowEdge(i)>=min:
            if (xmin>1e5): xmin = hist.GetXaxis().GetBinLowEdge(i)
            if hist.GetXaxis().GetBinUpEdge(i)<=max:
                xmax = hist.GetXaxis().GetBinUpEdge(i)
                Nbins += 1
    return (Nbins,xmin,xmax)

def GetConfidenceIntervals(func, fitRes, Nbins, xcenters, ci, cl = 0.68):
    npar = func.GetNumberFreeParameters()
    npar_real = func.GetNpar()
    fixed = 0
    if npar_real != npar:
        print "ERROR. Check other function."

    covmatr = fitRes.GetCovarianceMatrix()
    rho = fitRes.GetCorrelationMatrix()
    tStudent = ROOT.TMath.StudentQuantile(0.5 + cl/2, func.GetNDF())
    chindf = ROOT.TMath.Sqrt(func.GetChisquare()/func.GetNDF())

    grad = array('d',[0.]*npar_real)
    sum_vector = array('d',[0.]*npar)
    for ipoint in range(Nbins):
        ci_=0
        func.GradientPar(array('d',[xcenters[ipoint]]), grad)
        # multiply the covariance matrix by gradient
        for irow in range(npar):
            sum_vector[irow]=0;
            for icol in range(npar):
                igrad=0
                ifree=0
                if fixed:
                    print "ERROR. Check other function."
                else:
                    igrad = icol;
                sum_vector[irow] += covmatr[irow][icol]*grad[igrad]
        igrad = 0;
        for i_ in range(npar):
            igrad=0
            ifree=0
            if fixed:
                print "ERROR. Check other function."
            else:
                igrad = i_
            ci_ += grad[igrad]*sum_vector[i_]
        ci_ = ROOT.TMath.Sqrt(ci_)
        ci[ipoint] = ci_*tStudent*chindf

def ComputeHistWithCL(name, func, fitRes, hist, cl=0.68):
  # create a histogram for the fit region only
  Nbins,xmin,xmax = CountBinsMinMax(hist, func.GetXmin(), func.GetXmax())
  band = ROOT.TH1F("band"+name+str(cl),name+str(cl), Nbins, xmin, xmax)
  band_pull  = ROOT.TH1F("band_pull"+name+str(cl),name+str(cl), Nbins, xmin, xmax)
  pull = ROOT.TH1F("pull"+name+str(cl),name+str(cl), Nbins, xmin, xmax)

  # now get an array with the bin centers and compute the CLs
  xcenters = array('d',[0.]*Nbins)
  ci = array('d',[0.]*Nbins)
  for i in range(1,Nbins+1):
      xcenters[i-1] = band.GetXaxis().GetBinCenter(i)

  GetConfidenceIntervals(func, fitRes, Nbins, xcenters, ci, cl);

  for i in range(1,Nbins+1):
      x_ = band.GetXaxis().GetBinCenter(i)
      y_ = func.Eval(x_)
      if y_==0: y_ = 1e-07
      ibin = hist.GetXaxis().FindBin(x_)
      band.SetBinContent(i, y_)
      band.SetBinError(i, ci[i-1])
      band_pull.SetBinContent(i, 1)
      band_pull.SetBinError(i, ci[i-1]/y_)
      pull.SetBinContent(i, hist.GetBinContent(ibin)/y_)
      pull.SetBinError(i, hist.GetBinError(ibin)/y_)
      # print name, x_, pull.GetBinContent(i), band_pull.GetBinError(i)
  return band, band_pull,pull



def CrystalBall_Fit(x, par):
    mean = par[0]
    sigma = par[1]
    alpha = par[2]
    n = par[3]
    norm = par[4]
    std =(x[0]-mean)/sigma
    if (alpha < 0): std *= -1
    alpha = ROOT.TMath.Abs(alpha)
    result = 0
    if std >= -alpha:
        result = ROOT.TMath.Exp(-0.5*std*std)
    else:
        A = ROOT.TMath.Power(n/alpha, n)*ROOT.TMath.Exp(-0.5*alpha*alpha)
        B = n/alpha-alpha
        result = A/ROOT.TMath.Power(B-std, n)
    return norm*result


def ExpGaussExp_Fit(x, par):
    mean = par[0]
    sigma = par[1]
    kl = par[2]
    kh = par[3]
    norm = par[4]
    std =(x[0]-mean)/sigma
    result = 0
    if std < -kl:
        result = ROOT.TMath.Exp(kl*kl*0.5+kl*std)
    elif -kl <= std < kh:
        result = ROOT.TMath.Exp(-0.5*std*std)
    else:
        result = ROOT.TMath.Exp(kh*kh*0.5-kh*std)
    return norm*result
