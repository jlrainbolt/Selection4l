#include <tuple>

#include "TColor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMathText.h"





///////////////////////
//                   //
//    DECLARATION    //
//                   //
///////////////////////


//--- FACELIFT FUNCTIONS ---//

void Facelift(TCanvas*);
void Facelift(TH1*);
void Facelift(THStack*);    // MUST CALL AFTER DRAWING
void Facelift(TLegend*);    // MUST CALL AFTER DRAWING



//--- PLOT DIMENSIONS ---//

const UInt_t    lCanvasSize = 800;                  // Canvas width/height
const Float_t   lCanvasMargin = 0.12;               // Give us more room for titles

const Font_t    lHelveticaMediumR = 43;             // Font code 3 is expressed in pixels
const Float_t   lLarge = 36, lSmall = 24;           // Font sizes
const Float_t   lTitleOffsetY = 1.1;                // Fix for y axis offset

const Float_t   lLegendMargin = 0.04;               // Margin for legend from upper right corner



//--- COLORMAP INTERFACE ---//

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

// "Transparency" utilities
tuple<float, float, float> GetAlpha(const tuple<float, float, float>, const float);
float GetAlpha(const float, const float);



//--- LaTeX COMMANDS ---//

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

TString _BFZto(const TString chan)  {return "\\mathcal{B}(\\" + _Z + _to + chan + ")";}
TString _pT(const TString x)        {return "p_{\\mbox{T} " + x + "}";}
TString _pT(const unsigned n)       {return _pT(TString::Format("%u", n));}
TString _p4(const TString x)        {return "\\vec{p}_{\\mbox{T} " + x + "}";}
TString _p4(const unsigned n)       {return _p4(TString::Format("%u", n));}
TString _m(const TString x)         {return "m_{" + x + "}";}
TString _m(const unsigned n)        {return _m(TString::Format("%u", n));}
TString _eta(const TString x)       {return "\\eta_{" + x + "}";}
TString _eta(const unsigned n)      {return _eta(TString::Format("%u", n));}

const TString   _GeV    = "\\mbox{GeV}";
const TString   _rad    = "\\mbox{rad}";
const TString   _unit   = "\\mbox{unit}";
const TString   _units  = "\\mbox{units}";

TString _EventsPer(const TString x, const TString unit) {return "\\mbox{Events / }" + x + unit;}
TString _EventsPer(const TString unit)                  {return _EventsPer("", unit);}
/*
TString _EventsPer(const unsigned n, const TString unit) {
    return _EventsPer(TString::Format("%u", n), unit);
}
*/
TString _EventsPer(const float n, const TString unit) {
    return _EventsPer(TString::Format("%g", n), unit);
}

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



/////////////////////////
//                     //
//    MAIN FUNCTION    //
//                     //
/////////////////////////


void loadPlotSetup()
{


    //--- CREATE COLORS ---//

    int ColorIndex[7] = {lBlue, lOrange, lYellow, lPurple, lGreen, lLightBlue, lRed};
    tuple<float, float, float> Lines[7] = { std::make_tuple(0, 0.4470, 0.7410),
                                            std::make_tuple(0.8500, 0.3250, 0.0980),
                                            std::make_tuple(0.9290, 0.6940, 0.1250),
                                            std::make_tuple(0.4940, 0.1840, 0.5560),
                                            std::make_tuple(0.4660, 0.6740, 0.1880),
                                            std::make_tuple(0.3010, 0.7450, 0.9330),
                                            std::make_tuple(0.6350, 0.0780, 0.1840) };

    for (unsigned i = 0; i < 7; i++)
    {
        auto rgb = GetAlpha(Lines[i], lAlpha);
        float r, g, b;
        std::tie(r, g, b) = rgb;
        
        lines[i] = new TColor(ColorIndex[i], r, g, b);
    }



    //--- FORCE TMathText ---//




    cout << "Loaded colors, text commands, and facelift functions for plots" << endl;    
}





//////////////////////////
//                      //
//    PLOT FACELIFTS    //
//                      //
//////////////////////////



//--- TCanvas ---//

void Facelift(TCanvas *canvas)
{
    canvas->SetCanvasSize(lCanvasSize, lCanvasSize);
    canvas->SetMargin(lCanvasMargin, lCanvasMargin, lCanvasMargin, lCanvasMargin);
}



//--- TH1 ---//

void Facelift(TH1 *hist)
{
    hist->SetStats(0);
    hist->SetMinimum(0);

    hist->GetXaxis()->SetLabelFont(lHelveticaMediumR);
    hist->GetXaxis()->SetLabelSize(lSmall);
    hist->GetXaxis()->SetTitleFont(lHelveticaMediumR);
    hist->GetXaxis()->SetTitleSize(lLarge);

    hist->GetYaxis()->SetLabelFont(lHelveticaMediumR);
    hist->GetYaxis()->SetLabelSize(lSmall);
    hist->GetYaxis()->SetTitleFont(lHelveticaMediumR);
    hist->GetYaxis()->SetTitleSize(lLarge);
    hist->GetYaxis()->SetTitleOffset(lTitleOffsetY);
}



//--- THStack ---//

void Facelift(THStack *stack)
{
    // MUST CALL AFTER DRAWING

    stack->SetMinimum(0);

    stack->GetXaxis()->SetLabelFont(lHelveticaMediumR);
    stack->GetXaxis()->SetLabelSize(lSmall);
    stack->GetXaxis()->SetTitleFont(lHelveticaMediumR);
    stack->GetXaxis()->SetTitleSize(lLarge);

    stack->GetYaxis()->SetLabelFont(lHelveticaMediumR);
    stack->GetYaxis()->SetLabelSize(lSmall);
    stack->GetYaxis()->SetTitleFont(lHelveticaMediumR);
    stack->GetYaxis()->SetTitleSize(lLarge);
    stack->GetYaxis()->SetTitleOffset(lTitleOffsetY);

//  gPad->Modified();
}



//--- TLegend ---//

void Facelift(TLegend *legend)
{
    // MUST CALL AFTER DRAWING

    legend->SetTextFont(lHelveticaMediumR);
    legend->SetTextSize(lLarge);
/*
    float LeftPosition = 0.5,       LeftMargin = 2. * lCanvasMargin - lLegendMargin;
    float RightPosition = 1,        RightMargin = -lLegendMargin;
    float TopPosition = 1,          TopMargin = -lLegendMargin;
    float BottomPosition = TopPosition - 0.065 * (float) legend->GetNRows();
    float BottomMargin = 2. * lCanvasMargin - lLegendMargin;

    legend->SetX1NDC(LeftPosition + LeftMargin);
    legend->SetX2NDC(TopPosition + TopMargin);
    legend->SetY1NDC(BottomPosition - TopMargin);
    legend->SetY2NDC(TopPosition + TopMargin);

    gPad->Modified();
*/
}





///////////////////////
//                   //
//    COLOR UTILS    //
//                   //
///////////////////////


// Fake transparency by changing RGB
tuple<float, float, float> GetAlpha(const tuple<float, float, float> rgb, const float alpha)
{
    float r, g, b;
    std::tie(r, g, b) = rgb;
    return std::make_tuple(GetAlpha(r, alpha), GetAlpha(g, alpha), GetAlpha(b, alpha));
}


float GetAlpha(const float x, const float alpha)
{
    return x + (1 - x) * (1 - alpha);
}
