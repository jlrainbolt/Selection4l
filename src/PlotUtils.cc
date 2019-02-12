// STL
#include <tuple>

// ROOT
#include "TColor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"

// Custom
#include "PlotUtils.hh"

using namespace std;



//
//  FACELIFT
//


// FIXME reduce to TAxis function


// TCanvas
void Facelift(TCanvas *canvas)
{
    canvas->SetCanvasSize(lCanvasSize, lCanvasSize);
    canvas->SetMargin(lCanvasMargin, lCanvasMargin, lCanvasMargin, lCanvasMargin);
}

// TH1
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

// TAxis
void Facelift(TAxis *axis)
{
    axis->SetLabelFont(lHelveticaMediumR);
    axis->SetLabelSize(lSmall);
    axis->SetTitleFont(lHelveticaMediumR);
    axis->SetTitleSize(lLarge);
}

// THStack
void Facelift(THStack *stack)
{
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
}

// TLegend
void Facelift(TLegend *legend)
{
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



//
//  COLORS
//

// GetAlpha
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
