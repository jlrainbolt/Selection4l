#ifndef PLOTUTILS_HH
#define PLOTUTILS_HH


// ROOT
#include "TColor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"



//
//  PLOT DIMENSIONS
//

const UInt_t    lCanvasSize = 800;                  // Canvas width/height
const Float_t   lCanvasMargin = 0.12;               // Give us more room for titles

const Font_t    lHelveticaMediumR = 43;             // Font code 3 is expressed in pixels
const Float_t   lLarge = 36, lSmall = 24;           // Font sizes
const Float_t   lTitleOffsetY = 1.1;                // Fix for y axis offset

const Float_t   lLegendMargin = 0.04;               // Margin for legend from upper right corner



//
//  COLOR INTERFACE
//

TColor *lines[7];

const Int_t     lFreeColorIndex = 1179;             // First empty color index
const Int_t     lBlue       = lFreeColorIndex + 0;
const Int_t     lOrange     = lFreeColorIndex + 1;
const Int_t     lYellow     = lFreeColorIndex + 2;
const Int_t     lPurple     = lFreeColorIndex + 3;
const Int_t     lGreen      = lFreeColorIndex + 4;
const Int_t     lLightBlue  = lFreeColorIndex + 5;
const Int_t     lRed        = lFreeColorIndex + 6;

const Float_t   lAlpha      = 0.75;



//
// "LaTeX" COMMANDS
//

const TString   _l      = "\\ell";
const TString   _mu     = "\\mu";
const TString   _e      = "\\mbox{e}";

const TString   _to     = "\\to";
const TString   _plus   = "^{+}";
const TString   _minus  = "^{-}";

const TString   _ll     = _l + _plus + _l + _minus;
const TString   _mumu   = _mu + _plus + _mu + _minus;
const TString   _ee     = _e + _plus + _e + _minus;
const TString   _4l     = "4" + _l;
const TString   _4mu    = "4" + _mu;
const TString   _2mu2e  = "2" + _mu + "2" + _e;
const TString   _2e2mu  = "2" + _e + "2" + _mu;
const TString   _4e     = "4" + _e;

const TString   _Z      = "\\mbox{Z}";
const TString   _ZZ     = _Z + _Z;
const TString   _Z1     = _Z + "_{1}";
const TString   _Z2     = _Z + "_{2}";

const TString   _GeV    = "\\mbox{GeV}";
const TString   _rad    = "\\mbox{rad}";
const TString   _unit   = "\\mbox{unit}";
const TString   _units  = "\\mbox{units}";

const TString   _met    = "p_{\\mbox{T}}^{\\mbox{miss}}";
const TString   _iso    = "\\mbox{RelPFIso}";

const TString   _g      = "\\mbox{g}";
const TString   _q      = "\\mbox{q}";
const TString   _t      = "\\mbox{t}";
const TString   _tbar   = "\\bar{" + _t + "}";
const TString   _ttbar  = _t + _tbar;
const TString   _H      = "\\mbox{H}";
const TString   _W      = "\\mbox{W}";
const TString   _V      = "\\mbox{V}";
const TString   _nu     = "\\nu";
const TString   _jets   = "+ \\mbox{jets}";
const TString   _and    = "\\mbox{,}";

const TString   _sp     = "\\,";


TString _BFZto( const TString   chan)   {return "\\mathcal{B}(\\" + _Z + _to + chan + ")";}
TString _pT(    const TString   x)      {return "p_{\\mbox{T} " + x + "}";}
TString _pT(    const unsigned  n)      {return _pT(TString::Format("%u", n));}
TString _p4(    const TString   x)      {return "\\vec{p}_{\\mbox{T} " + x + "}";}
TString _p4(    const unsigned  n)      {return _p4(TString::Format("%u", n));}
TString _m(     const TString   x)      {return "m_{" + x + "}";}
TString _m(     const unsigned  n)      {return _m(TString::Format("%u", n));}
TString _eta(   const TString   x)      {return "\\eta_{" + x + "}";}
TString _eta(   const unsigned  n)      {return _eta(TString::Format("%u", n));}

/*
TString _EventsPer(const TString x, const TString unit)
    {return "\\mbox{Events / }" + x + "\\mbox{ }" + unit;}
TString _EventsPer(const TString unit)                  {return _EventsPer("", unit);}
TString _EventsPer(const unsigned n, const TString unit) {
    return _EventsPer(TString::Format("%u", n), unit);
}
TString _EventsPer(const float n, const TString unit) {
    return _EventsPer(TString::Format("%g", n), unit);
}
*/


//
//  FUNCTIONS
//

// Facelift
void Facelift(TCanvas*);
void Facelift(TH1*);
void Facelift(THStack*);    // Must call after drawing
void Facelift(TLegend*);    // Must call after drawing

// "Transparency" utilities
tuple<float, float, float> GetAlpha(const tuple<float, float, float>, const float);
float GetAlpha(const float, const float);



#endif  // PLOTUTILS_HH