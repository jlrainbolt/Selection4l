#include <tuple>

#include "TColor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMathText.h"



//--- HELPERS ---//

void Facelift(TCanvas*);
void Facelift(TH1*);
void Facelift(THStack*);    // MUST CALL AFTER DRAWING
void Facelift(TLegend*);    // MUST CALL AFTER DRAWING

tuple<float, float, float> GetAlpha(const tuple<float, float, float>, const float);
float GetAlpha(const float, const float);



//--- GLOBAL VARIABLES ---//

// Histogram sizing
const UInt_t    lCanvasSize = 800;                  // Canvas width/height
const Float_t   lCanvasMargin = 0.12;               // Give us more room for titles

const Font_t    lHelveticaMediumR = 43;             // Font code 3 is expressed in pixels
const Float_t   lLarge = 36, lSmall = 24;           // Font sizes
const Float_t   lTitleOffsetY = 1.1;                // Fix for y axis offset

const Float_t   lLegendMargin = 0.04;               // Margin for legend from upper right corner


// Colormap
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



// LaTeX strings
const TString   _l      = "\\ell";
const TString   _mu     = "\\mu";
const TString   _e      = "\\mbox{e}";
const TString   _plus   = "^{+}";
const TString   _minus  = "^{-}";

const TString   _ll     = _l + _plus + _l + _minus;
const TString   _mumu   = _mu + _plus + _mu + _minus;
const TString   _ee     = _e + _plus + _e + _minus;
const TString   _4l     = "4" + _ell;
const TString   _4mu    = "4" + _mu;
const TString   _2mu2e  = "2" + _mu + "2" + _e;
const TString   _2e2mu  = "2" + _e + "2" + _mu;
const TString   _4e     = "4" + _e;

const TString   _Z      = "\\mbox{Z}";
const TString   _Z      = _Z + _Z;
const TString   _Z1     = _Z + "_{1}";
const TString   _Z2     = _Z + "_{2}";

TString _BFZto(const TString chan)  {return "\\mathcal{B}(\\" + Z + " \\to " + decay + ")";}
TString _pT(const TString x)        {return "p_{\\mbox{T} " + x + "}";}
TString _pT(const unsigned n)       {return _pT(TString::Format("%u", n));}
TString _m(const TString x)         {return "m_{" + x + "}";}
TString _m(const unsigned n)        {return _m(TString::Format("%u", n));}

const TString   _GeV    = "\\mbox{ GeV}";
const TString   _rad    = "\\mbox{ rad}";
const TString   _unit   = "\\mbox{ unit}";
const TString   _units  = "\\mbox{ units}";

TString _EventsPer(const TString x, const TString unit) {return "\\mbox{Events / }" + x + unit;}
TString _EventsPer(const float n, const TString unit) {
    return TString::Format("\\mbox{Events / } %g" + unit, n);   }



//--- MAIN FUNCTION ---//

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


    // Force TMathText rendering
//  TString "\\star";


    cout << "Loaded custom colors and facelift functions for plots" << endl;    
}




void Facelift(TCanvas *canvas)
{
    canvas->SetCanvasSize(lCanvasSize, lCanvasSize);
    canvas->SetMargin(lCanvasMargin, lCanvasMargin, lCanvasMargin, lCanvasMargin);
}


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

    gPad->Modified();
}


void Facelift(TLegend *legend)
{
    // MUST CALL AFTER DRAWING

    legend->SetTextFont(lHelveticaMediumR);
    legend->SetTextSize(lLarge);

    float LeftPosition = 0.5,       LeftMargin = 2. * lCanvasMargin - lLegendMargin;
    float RightPosition = 1,        RightMargin = -lLegendMargin;
    float TopPosition = 1,          TopMargin = -lLegendMargin;
    float BottomPosition = TopPosition - 0.065 * (float) legend->GetNRows();
    float BottomMargin = 2. * lCanvasMargin - lLegendMargin;

    legend->SetX1NDC(LeftPosition + LeftMargin);
    legend->SetX2NDC(TopPosition + TopMargin);
//  legend->SetY1NDC(BottomPosition + TopMargin);
    legend->SetY1NDC(0.5 + BottomMargin);
    legend->SetY2NDC(TopPosition + TopMargin);

    gPad->Modified();
}


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
