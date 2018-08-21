#include "PlotUtils.h"


void PlotUtils()
{
    // Create colors from Lines colormap
    TColor *lines[7];
    for (unsigned i = 0; i < 7; i++)
    {
        auto rgb = GetAlpha(Lines[i], lAlpha);
        float r, g, b;
        std::tie(r, g, b) = rgb;
        
        lines[i] = new TColor(lFreeColorIndex + (int) i, r, g, b);
    }


    // Force TMathText rendering
//  TString "\\star";
    
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


void Facelift(TLegend *legend)
{
    legend->SetTextFont(lHelveticaMediumR);
    legend->SetTextSize(lLarge);
   

    float BottomLeftMargin = 2. * lCanvasMargin - lLegendMargin, TopRightMargin = lLegendMargin;
    float LeftPosition = 0.5, TopRightPosition = 1;
    float BottomPosition = 0.125 * (float) legend->GetNRows();

    legend->SetX1(LeftPosition + BottomLeftMargin);
    legend->SetY1(BottomPosition + BottomLeftMargin);
    legend->SetX2(TopRightPosition - TopRightMargin);
    legend->SetY2(TopRightPosition - TopRightMargin);
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
