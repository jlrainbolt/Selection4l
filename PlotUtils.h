#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include <tuple>

#include "TColor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMathText.h"


// Functions
void Facelift(TCanvas*);
void Facelift(TH1*);
void Facelift(THStack*);
void Facelift(TLegend*);

tuple<float, float, float> GetAlpha(const tuple<float, float, float>, const float);
float GetAlpha(const float, const float);


// "Data members"
UInt_t  lCanvasSize = 800;                  // Canvas width/height
Float_t lCanvasMargin = 0.12;               // Give us more room for titles

Font_t  lHelveticaMediumR = 43;             // Font code 3 is expressed in pixels
Float_t lLarge = 36, lSmall = 24;           // Font sizes
Float_t lTitleOffsetY = 1.1;                // Fix for y axis offset

Float_t lLegendMargin = 0.04;               // Margin for legend from upper right corner


// Colormap
tuple<Float_t, Float_t, Float_t> Lines[7] = {   std::make_tuple(0, 0.4470, 0.7410),
                                                std::make_tuple(0.8500, 0.3250, 0.0980),
                                                std::make_tuple(0.9290, 0.6940, 0.1250),
                                                std::make_tuple(0.4940, 0.1840, 0.5560),
                                                std::make_tuple(0.4660, 0.6740, 0.1880),
                                                std::make_tuple(0.3010, 0.7450, 0.9330),
                                                std::make_tuple(0.6350, 0.0780, 0.1840) };

Int_t   lFreeColorIndex = 1179;             // First empty color index
Int_t   lBlue       = lFreeColorIndex + 0;
Int_t   lOrange     = lFreeColorIndex + 1;
Int_t   lYellow     = lFreeColorIndex + 2;
Int_t   lPurple     = lFreeColorIndex + 4;
Int_t   lGreen      = lFreeColorIndex + 5;
Int_t   lLightBlue  = lFreeColorIndex + 6;
Int_t   lRed        = lFreeColorIndex + 7;

Float_t lAlpha      = 0.75;

#endif
