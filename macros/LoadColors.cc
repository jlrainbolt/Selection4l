// ROOT
#include "TColor.h"
#include "TMathText.h"

// Custom
#include "PlotUtils.hh"



void LoadColors()
{
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

    cout << "Loaded colors from PlotUtils.hh" << endl; 
}
