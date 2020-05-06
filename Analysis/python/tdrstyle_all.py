import ROOT as rt

#######################################
# New CMS Style from 2014             #
# https://ghm.web.cern.ch/ghm/plots/  #
# Merged all macros into one          #
#######################################

#############
# tdrstyle  #
#############

def tdrGrid( gridOn):
  tdrStyle.SetPadGridX(gridOn)
  tdrStyle.SetPadGridY(gridOn)

#fixOverlay: Redraws the axis
def fixOverlay():
  gPad.RedrawAxis()

def setTDRStyle():
  tdrStyle =  rt.TStyle("tdrStyle","Style for P-TDR")
  #for the canvas:
  tdrStyle.SetCanvasBorderMode(0)
  tdrStyle.SetCanvasColor(rt.kWhite)
  tdrStyle.SetCanvasDefH(600) #Height of canvas
  tdrStyle.SetCanvasDefW(600) #Width of canvas
  tdrStyle.SetCanvasDefX(0)   #POsition on screen
  tdrStyle.SetCanvasDefY(0)
  tdrStyle.SetPadBorderMode(0)
  #tdrStyle.SetPadBorderSize(Width_t size = 1)
  tdrStyle.SetPadColor(rt.kWhite)
  tdrStyle.SetPadGridX(False)
  tdrStyle.SetPadGridY(False)
  tdrStyle.SetGridColor(0)
  tdrStyle.SetGridStyle(3)
  tdrStyle.SetGridWidth(1)
  #For the frame:
  tdrStyle.SetFrameBorderMode(0)
  tdrStyle.SetFrameBorderSize(1)
  tdrStyle.SetFrameFillColor(0)
  tdrStyle.SetFrameFillStyle(0)
  tdrStyle.SetFrameLineColor(1)
  tdrStyle.SetFrameLineStyle(1)
  tdrStyle.SetFrameLineWidth(1)
  #For the histo:
  #tdrStyle.SetHistFillColor(1)
  #tdrStyle.SetHistFillStyle(0)
  tdrStyle.SetHistLineColor(1)
  tdrStyle.SetHistLineStyle(0)
  tdrStyle.SetHistLineWidth(1)
  #tdrStyle.SetLegoInnerR(Float_t rad = 0.5)
  #tdrStyle.SetNumberContours(Int_t number = 20)
  tdrStyle.SetEndErrorSize(2)
  #tdrStyle.SetErrorMarker(20)
  #tdrStyle.SetErrorX(0.)
  tdrStyle.SetMarkerStyle(20)
  #For the fit/function:
  tdrStyle.SetOptFit(1)
  tdrStyle.SetFitFormat("5.4g")
  tdrStyle.SetFuncColor(2)
  tdrStyle.SetFuncStyle(1)
  tdrStyle.SetFuncWidth(1)
  #For the date:
  tdrStyle.SetOptDate(0)
  # tdrStyle.SetDateX(Float_t x = 0.01)
  # tdrStyle.SetDateY(Float_t y = 0.01)
  # For the statistics box:
  tdrStyle.SetOptFile(0)
  tdrStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
  tdrStyle.SetStatColor(rt.kWhite)
  tdrStyle.SetStatFont(42)
  tdrStyle.SetStatFontSize(0.025)
  tdrStyle.SetStatTextColor(1)
  tdrStyle.SetStatFormat("6.4g")
  tdrStyle.SetStatBorderSize(1)
  tdrStyle.SetStatH(0.1)
  tdrStyle.SetStatW(0.15)
  # tdrStyle.SetStatStyle(Style_t style = 1001)
  # tdrStyle.SetStatX(Float_t x = 0)
  # tdrStyle.SetStatY(Float_t y = 0)
  # Margins:
  tdrStyle.SetPadTopMargin(0.05)
  tdrStyle.SetPadBottomMargin(0.13)
  tdrStyle.SetPadLeftMargin(0.16)
  tdrStyle.SetPadRightMargin(0.02)
  # For the Global title:
  tdrStyle.SetOptTitle(0)
  tdrStyle.SetTitleFont(42)
  tdrStyle.SetTitleColor(1)
  tdrStyle.SetTitleTextColor(1)
  tdrStyle.SetTitleFillColor(10)
  tdrStyle.SetTitleFontSize(0.05)
  # tdrStyle.SetTitleH(0) # Set the height of the title box
  # tdrStyle.SetTitleW(0) # Set the width of the title box
  # tdrStyle.SetTitleX(0) # Set the position of the title box
  # tdrStyle.SetTitleY(0.985) # Set the position of the title box
  # tdrStyle.SetTitleStyle(Style_t style = 1001)
  # tdrStyle.SetTitleBorderSize(2)
  # For the axis titles:
  tdrStyle.SetTitleColor(1, "XYZ")
  tdrStyle.SetTitleFont(42, "XYZ")
  tdrStyle.SetTitleSize(0.06, "XYZ")
  # tdrStyle.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
  # tdrStyle.SetTitleYSize(Float_t size = 0.02)
  tdrStyle.SetTitleXOffset(0.9)
  tdrStyle.SetTitleYOffset(1.25)
  # tdrStyle.SetTitleOffset(1.1, "Y") # Another way to set the Offset
  # For the axis labels:
  tdrStyle.SetLabelColor(1, "XYZ")
  tdrStyle.SetLabelFont(42, "XYZ")
  tdrStyle.SetLabelOffset(0.007, "XYZ")
  tdrStyle.SetLabelSize(0.05, "XYZ")
  # For the axis:
  tdrStyle.SetAxisColor(1, "XYZ")
  tdrStyle.SetStripDecimals(True)
  tdrStyle.SetTickLength(0.03, "XYZ")
  tdrStyle.SetNdivisions(510, "XYZ")
  tdrStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
  tdrStyle.SetPadTickY(1)
  # Change for log plots:
  tdrStyle.SetOptLogx(0)
  tdrStyle.SetOptLogy(0)
  tdrStyle.SetOptLogz(0)
  # Postscript options:
  tdrStyle.SetPaperSize(20.,20.)
  # tdrStyle.SetLineScalePS(Float_t scale = 3)
  # tdrStyle.SetLineStyleString(Int_t i, const char* text)
  # tdrStyle.SetHeaderPS(const char* header)
  # tdrStyle.SetTitlePS(const char* pstitle)
  # tdrStyle.SetBarOffset(Float_t baroff = 0.5)
  # tdrStyle.SetBarWidth(Float_t barwidth = 0.5)
  # tdrStyle.SetPaintTextFormat(const char* format = "g")
  # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
  # tdrStyle.SetTimeOffset(Double_t toffset)
  # tdrStyle.SetHistMinimumZero(kTRUE)
  tdrStyle.SetHatchesLineWidth(5)
  tdrStyle.SetHatchesSpacing(0.05)
  tdrStyle.cd()

###################
# end setTDRStyle #
###################
#################
# end tdrstyle  #
#################


#############
# CMS_lumi  #
#############

#
# Global variables
#

cmsText     = "CMS"
cmsTextFont   = 61

writeExtraText = True
extraText   = "Preliminary"
extraText2  = "" # For Simulation Preliminary on two lines
extraTextFont = 52

# text sizes and text offsets with respect to the top frame in unit of the top margin size
lumiTextSize     = 0.6
lumiTextOffset   = 0.2
cmsTextSize      = 0.75
cmsTextOffset    = 0.1

relPosX    = 0.045
relPosY    = 0.035
relExtraDY = 1.2

# ratio of "CMS" and extra text size
extraOverCmsTextSize  = 0.76

lumi_13TeV = "41.53 fb^{-1}"
lumi_8TeV  = "19.7 fb^{-1}"
lumi_7TeV  = "5.1 fb^{-1}"
lumi_sqrtS = ""

drawLogo = False

def CMS_lumi(pad, iPeriod=4, iPosX=11):
  outOfFrame    = False
  if(iPosX/10==0 ): outOfFrame = True
  alignY_=3
  alignX_=2
  if( iPosX/10==0 ): alignX_=1
  if( iPosX==0    ): alignX_=1
  if( iPosX==0    ): alignY_=1
  if( iPosX/10==1 ): alignX_=1
  if( iPosX/10==2 ): alignX_=2
  if( iPosX/10==3 ): alignX_=3
  # if( iPosX==0    ): relPosX = pad.GetLeftMargin()
  align_ = 10*alignX_ + alignY_
  H = pad.GetWh()
  W = pad.GetWw()
  l = pad.GetLeftMargin()
  t = pad.GetTopMargin()
  r = pad.GetRightMargin()
  b = pad.GetBottomMargin()
  e = 0.025
  pad.cd()
  lumiText = ""
  if( iPeriod==1 ):
    lumiText += lumi_7TeV
    lumiText += " (7 TeV)"
  elif ( iPeriod==2 ):
    lumiText += lumi_8TeV
    lumiText += " (8 TeV)"
  elif( iPeriod==3 ):
    lumiText = lumi_8TeV
    lumiText += " (8 TeV)"
    lumiText += " + "
    lumiText += lumi_7TeV
    lumiText += " (7 TeV)"
  elif ( iPeriod==4 ):
    lumiText += lumi_13TeV
    lumiText += " (13 TeV)"
  elif ( iPeriod==7 ):
    if( outOfFrame ):lumiText += "#scale[0.85]{"
    lumiText += lumi_13TeV
    lumiText += " (13 TeV)"
    lumiText += " + "
    lumiText += lumi_8TeV
    lumiText += " (8 TeV)"
    lumiText += " + "
    lumiText += lumi_7TeV
    lumiText += " (7 TeV)"
    if( outOfFrame): lumiText += "}"
  elif ( iPeriod==12 ):
    lumiText += "8 TeV"
  elif ( iPeriod==0 ):
    lumiText += lumi_sqrtS
  # print lumiText
  latex = rt.TLatex()
  latex.SetNDC()
  latex.SetTextAngle(0)
  latex.SetTextColor(rt.kBlack)
  extraTextSize = extraOverCmsTextSize*cmsTextSize
  latex.SetTextFont(42)
  latex.SetTextAlign(31)
  latex.SetTextSize(lumiTextSize*t)
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)
  if( outOfFrame ):
    latex.SetTextFont(cmsTextFont)
    latex.SetTextAlign(11)
    latex.SetTextSize(cmsTextSize*t)
    latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText)
  pad.cd()
  posX_ = 0
  if( iPosX%10<=1 ):   posX_ = l + relPosX*(1-l-r)
  elif( iPosX%10==2 ): posX_ = l + 0.5*(1-l-r)
  elif( iPosX%10==3 ): posX_ = 1-r - relPosX*(1-l-r)
  posY_ = 1-t - relPosY*(1-t-b)
  if( not outOfFrame ):
    if( drawLogo ):
      posX_ =   l + 0.045*(1-l-r)*W/H
      posY_ = 1-t - 0.045*(1-t-b)
      xl_0 = posX_
      yl_0 = posY_ - 0.15
      xl_1 = posX_ + 0.15*H/W
      yl_1 = posY_
      CMS_logo = rt.TASImage("CMS-BW-label.png")
      pad_logo =  rt.TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 )
      pad_logo.Draw()
      pad_logo.cd()
      CMS_logo.Draw("X")
      pad_logo.Modified()
      pad.cd()
    else:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextSize(cmsTextSize*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, cmsText)
        if( writeExtraText ) :
            latex.SetTextFont(extraTextFont)
            latex.SetTextAlign(align_)
            latex.SetTextSize(extraTextSize*t)
            latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText)
            # For Simulation Preliminary
            if (extraText2!=""): latex.DrawLatex(posX_, posY_-relExtraDY*cmsTextSize*t - relExtraDY*extraTextSize*t, extraText2)
  elif( writeExtraText ):
    if( iPosX==0):
      posX_ =   l +  relPosX*(1-l-r)
      posY_ =   1-t+lumiTextOffset*t
    latex.SetTextFont(extraTextFont)
    latex.SetTextSize(extraTextSize*t)
    latex.SetTextAlign(align_)
    latex.DrawLatex(posX_, posY_, extraText)
  pad.Update()

#################
# end CMS_lumi  #
#################


#################
# General Macro #
#################

kSquare = True
kRectangular = False

# Give the macro an empty histogram for h.Draw("AXIS")
# Create h after calling setTDRStyle to get all the settings right
def tdrCanvas(canvName, x_min, x_max, y_min, y_max, nameXaxis, nameYaxis, square=kRectangular, iPeriod=4, iPos=11, extraText_ = "Preliminary"):
  # setTDRStyle()
  # writeExtraText = true       # if extra text
  extraText  = extraText_  # default extra text is "Preliminary"
  # lumi_8TeV  = "19.5 fb^{-1}" # default is "19.7 fb^{-1}"
  # lumi_7TeV  = "5.0 fb^{-1}"  # default is "5.1 fb^{-1}"
  # int iPeriod = 3    # 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
  # second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  # iPos=11 : top-left, left-aligned
  # iPos=33 : top-right, right-aligned
  # iPos=22 : center, centered
  # iPos=0  : out of frame (in exceptional cases)
  # mode generally :
  #   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
  #  if( iPos==0 ) relPosX = 0.12
  W = 600 if square else 800
  H = 600 if square else 600
  # Simple example of macro: plot with CMS name and lumi text
  #  (this script does not pretend to work in all configurations)
  # iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
  # For instance:
  #               iPeriod = 3 means: 7 TeV + 8 TeV
  #               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV
  # Initiated by: Gautier Hamel de Monchenault (Saclay)
  #
  W_ref = 600 if square else 800
  H_ref = 600 if square else 600
  # references for T, B, L, R
  T = 0.07*H_ref if square else 0.08*H_ref
  B = 0.13*H_ref if square else 0.12*H_ref
  L = 0.15*H_ref if square else 0.12*H_ref
  R = 0.05*H_ref if square else 0.04*H_ref
  #
  canv = rt.TCanvas(canvName,canvName,50,50,W,H)
  canv.SetFillColor(0)
  canv.SetBorderMode(0)
  canv.SetFrameFillStyle(0)
  canv.SetFrameBorderMode(0)
  canv.SetLeftMargin( L/W )
  canv.SetRightMargin( R/W )
  canv.SetTopMargin( T/H )
  canv.SetBottomMargin( B/H )
  # FOR JEC plots, prefer to keep ticks on both sides
  #canv.SetTickx(0)
  #canv.SetTicky(0)
  #
  # assert(h)
  h = canv.DrawFrame(x_min,y_min,x_max,y_max)
  h.GetYaxis().SetTitleOffset(1.25 if square else 1.0)
  h.GetXaxis().SetTitleOffset(1.0  if square else 0.9)
  h.GetXaxis().SetTitle(nameXaxis)
  h.GetYaxis().SetTitle(nameYaxis)
  h.Draw("AXIS")
  #
  # writing the lumi information and the CMS "logo"
  # CMS_lumi( canv, iPeriod, iPos )
  #
  canv.Update()
  canv.RedrawAxis()
  canv.GetFrame().Draw()
  #
  return canv

#################
# end tdrCanvas #
#################




def setNewTitle(name = "new title"):
  gStyle.SetOptTitle(0)
  title = rt.TPaveLabel(.11,.95,.35,.99,name,"brNDC")
  title.Draw()

def tdrDraw(h, opt, marker=rt.kFullCircle, mcolor=rt.kBlack, lstyle=rt.kSolid, lcolor=-1, fstyle=1001, fcolor=rt.kYellow+1):
  h.SetMarkerStyle(marker)
  h.SetMarkerColor(mcolor)
  h.SetLineStyle(lstyle)
  h.SetLineColor(mcolor if lcolor==-1 else lcolor)
  h.SetFillStyle(fstyle)
  h.SetFillColor(fcolor)
  h.Draw(opt+"SAME")

def tdrLeg(x1, y1, x2, y2, textSize=0.045, textFont=42, textColor=rt.kBlack):
  leg = rt.TLegend(x1, y1, x2, y2, "", "brNDC")
  leg.SetFillStyle(rt.kNone)
  leg.SetBorderSize(0)
  leg.SetTextSize(textSize)
  leg.SetTextFont(textFont)
  leg.SetTextColor(textColor)
  leg.Draw()
  return leg

def tdrHeader(leg, legTitle, textAlign=12, textSize=0.04, textFont=42, textColor=rt.kBlack, isToRemove = True):
  header = rt.TLegendEntry( 0, legTitle, "h" )
  header.SetTextFont(textFont)
  header.SetTextSize(textSize)
  header.SetTextAlign(textAlign)
  header.SetTextColor(textColor)
  if isToRemove:
    leg.SetHeader(legTitle,"C")
    leg.GetListOfPrimitives().Remove(leg.GetListOfPrimitives().At(0))
    leg.GetListOfPrimitives().AddAt(header,0)
  else:
    leg.GetListOfPrimitives().AddLast(header)
