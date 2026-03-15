"""
AnthroPy: A Python Library for Anthropometric and Physiological Modelling.

A comprehensive toolkit for anthropometric, body composition, and metabolic profiling. This module implements canonical equations and indices from established literature, enabling the calculation of a wide range of metrics based on user-provided measurements.

Date created: 10th March 2026


Design
------
- The central entry point is `profiler(measurements)`, which accepts a dictionary of measurement values (weight, height, circumferences, bone breadths, skinfolds, age, sex, etc.).
- Each metric is computed only if its required parameters are present.
- Results are returned as a dictionary with metric names as keys and numeric values (or None if unavailable).
- Average values are computed within categories when multiple methods are available (e.g., BF-Average, BMR-Average, BSA-Average, LBM-Average, IBW-Average).


Usage
-----
>>> from anthropy import profiler
>>> measurements = {"mWeight": 75, "mHeight": 173, "aAge": 36, "aGender": "male", "cWaist": 85, "cHip": 95}
>>> results = profiler(measurements)
>>> results["AT_BMI"]
25.06


Notes
-----
- All formulas are implemented according to canonical references in anthropometry, physiology, and nutrition science.
- Units: weight in kilograms, height and circumferences in centimeters, skinfolds in millimeters, age in years.
- Somatotype calculations require corrected girths (subtracting skinfold thickness ÷ 10).


The following measures are implemented:
1. Attributes (AT)
    (a) AT_Siri: Convert body density to body fat percentage using the Siri equation.
    (b) AT_Brozek: Convert body density to body fat percentage using the Brozek equation.
    (c) AT_BMI: Calculate Body Mass Index (BMI).
    (d) AT_PI: Calculate Ponderal Index (PI).
    (e) AT_TCR: Thigh-to-Calf Ratio
    (f) AT_FWR: Forearm-to-Wrist Ratio
    (g) AT_CHtR: Calf-to-Height Ratio
    (h) AT_NHtR: Neck-to-Height Ratio

2: Adiposity Index (AI)
    (a) AI_WHR: Waist-to-Hip Ratio
    (b) AI_WHtR: Waist-to-Height Ratio
    (c) AI_WTR: Waist-to-Thigh Ratio
    (d) AI_ABSI: A Body Shape Index
    (e) AI_ConicityIndex: Conicity Index
    (f) AI_WCR: Waist-to-Calf Ratio
    (g) AI_WHtPI: Waist-to-Height Power Index
    (h) AI_BRI: Body Roundness Index

3. Percent Body Fat by Circumferences / Girth Measurements (BFc)
    (a) BFc_USNavy: US Navy method.
    (b) BFc_YMCA: YMCA method.
    (c) BFc_mYMCA: Modified YMCA method.
    (d) BFc_CovertBailey: Covert Bailey method.
    (e) BFc_BehnkeWilmore: Behnke-Wilmore circumference method.
    (f) BFc_RFM: Relative Fat Mass.
    (g) BFc_BAI: Body Adiposity Index.

4. Percent Body Fat by Skinfold Measurements (BFs)
    (a) BFs_JacksonPollock3: Jackson-Pollock 3-site method.
    (b) BFs_JacksonPollock4: Jackson-Pollock 4-site method.
    (c) BFs_JacksonPollock7: Jackson-Pollock 7-site method.
    (d) BFs_BehnkeWilmore: Behnke-Wilmore skinfold method.
    (e) BFs_DurninWomersley: Durnin–Womersley 4-site skinfold method.
    (f) BFs_Sloan: Sloan equation (1967)
    (g) BFs_Parrillo: Parrillo 9-site formula

5. Percent Body Fat by Mixed Measurements (BFx)
    (a) BFx_BMI: Estimate from BMI using Deurenberg's equation.
    (b) BFx_JacksonPollock3Girth: Jackson-Pollock 3-site + girth method.
    (c) BFx_Henry2018: Henry et al. (2018) method for Asian Chinese.

6. Basal Metabolic Rate (BMR)
    (a) BMR_MifflinStJeor: Mifflin-St Jeor equation.
    (b) BMR_HarrisBenedict: Harris-Benedict equation.
    (c) BMR_Kleiber: Kleiber's law.
    (d) BMR_RozaShizgal: Roza-Shizgal revision of the Harris-Benedict equation.
    (e) BMR_Schofield: Schofield equation.
    (f) BMR_Cunningham: Cunningham (1980) equation.
    (g) BMR_KatchMcArdle: Katch and McArdle (1983) equation.
    (h) BMR_FAO2004: FAO/WHO/UNU (2004) equations.

7: Body Surface Area (BSA)
    (a) BSA_DuBois: Du Bois & Du Bois (1916) formula.
    (b) BSA_Mosteller: Mosteller (1987) formula.
    (c) BSA_Haycock: Haycock et al. (1978) formula.
    (d) BSA_GehanGeorge: Gehan & George (1970) formula.
    (e) BSA_Fujimoto: Fujimoto formula

8: Lean Body Mass (LBM)
    (a) LBM_Boer: Boer (1984) formula.
    (b) LBM_James: James (1977) formula.
    (c) LBM_Hume: Hume (1966) formula.

9: Ideal Body Weight (IBW) 
    (a) IBW_Robinson: Robinson (1983) formula.
    (b) IBW_Miller: Miller et al. (1983) formula.
    (c) IBW_Hamwi: Hamwi (1964) formula.
    (d) IBW_Devine: Devine (1974) formula.
    (e) IBW_Lemmens: Lemmens et al. (2005) formula.
    (f) IBW_Peterson: Peterson et al. (2016) formula.
    (g) IBW_KeysBrozek: Keys & Brozek (1953) formula.

10: Sum of Skinfolds (SSF) defined for
    (a) Sum2 (sTricep, sSubscapular)
    (b) Sum3_men (sChest, sAbdominal, sThigh), equivalent to Jackson–Pollock for men
    (c) Sum3_women (sTricep, sSuprailiac, sThigh), equivalent to Jackson–Pollock for women
    (d) Sum4 (sTricep, sBicep, sSubscapular, sSuprailiac), equivalent to Durnin–Womersley
    (e) Sum5 (sTricep, sSubscapular, sSuprailiac, sAbdominal, sThigh)
    (f) Sum7 (sChest, sMidaxillary, sTricep, sSubscapular, sAbdominal, sSuprailiac, sThigh), equivalent to Jackson–Pollock 7-sites
    (g) Sum9 (sChest, sAbdominal, sThigh, sTricep, sSubscapular, sSuprailiac, sLowerback, sCalf, sBicep), equivalent to Parrillo 9-sites
    (h) All (sAbdominal, sBicep, sCalf, sChest, sLowerback, sMidaxillary, sSubscapular, sSuprailiac, sThigh, sTricep)
    (i) Any other combination by defining a list

11. Body Frame Classification (FC)
    (a) FC_Metropolitan: Classification using Metropolitan Life (1983)

12. Somatotype (ST)
    (a) ST_HeathCarter: Heath and Carter method.


License: GNU General Public License version 3 for academic or not-for-profit use only.

Bactome package is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import math


def _round(value: float, dp: int | None):
    return round(value, dp) if dp is not None else value


def AT_Siri(aDensity: float, dp: int = 2) -> float:
    """
    Convert body density to body fat percentage using the Siri equation.

    Reference: Siri, W.E. (1961). Body composition from fluid spaces and density: analysis of methods. In: Techniques for Measuring Body Composition. Washington, DC: National Academy of Sciences.

    Parameters:
        aDensity (float): Body density in g/cm³.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    if aDensity <= 0:
        raise ValueError("aDensity must be a positive number.")
    bf = (495 / aDensity) - 450
    return _round(bf, dp)


def AT_Brozek(aDensity: float,dp: int = 2) -> float:
    """
    Convert body density to body fat percentage using the Brozek equation.

    Reference: Brozek, J., Grande, F., Anderson, J.T. and Keys, A. (1963) Densitometric Analysis of Body Composition: Revision of Some Quantitative Assumption. Annals of the New York Academy of Sciences, 110, 113-140. https://doi.org/10.1111/j.1749-6632.1963.tb17079.x 

    Parameters:
        aDensity (float): Body density in g/cm³.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    if aDensity <= 0:
        raise ValueError("aDensity must be a positive number.")
    bf = ((4.57 / aDensity) - 4.142) * 100
    return _round(bf, dp)


def AT_BMI(mWeight: float, mHeight: float,dp: int = 2) -> float:
    """
    Calculate Body Mass Index (BMI).

    Reference: Keys, A., Fidanza, F., Karvonen, M. J., Kimura, N., & Taylor, H. L. (1972). Indices of relative weight and obesity. Journal of Chronic Diseases, 25(6), 329-343.

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: BMI value (kg/m²).
    """
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mHeight must be positive numbers.")
    # Convert height to meters
    mHeight_m = mHeight / 100
    bmi = mWeight / (mHeight_m * mHeight_m)
    return _round(bmi, dp)


def AT_PI(mWeight: float, mHeight: float,dp: int = 2) -> float:
    """
    Calculate Ponderal Index (PI).

    Reference: Davies, D. P. (1980). Size at birth and growth in the first year of life of babies who are overweight and underweight at birth. Proceedings of the Nutrition Society. 39 (1): 25–33. doi:10.1079/PNS19800005

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: BMI value (kg/m²).
    """
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mHeight must be positive numbers.")
    # Convert height to meters
    mHeight_m = mHeight / 100
    ponderal = mWeight / (mHeight_m * mHeight_m * mHeight_m)
    return _round(ponderal, dp)


def AT_TCR(cThigh: float, cCalf: float, dp: int = 4) -> float:
    """
    Calculate Thigh-to-Calf Ratio.

    Parameters:
        cThigh (float): Thigh circumference in cm.
        cCalf (float): Calf circumference in cm.
        dp (int): Number of decimal points. Default = 4.

    Returns:
        float: Thigh-to-Calf Ratio.
    """
    if cThigh <= 0 or cCalf <= 0:
        raise ValueError("cThigh and cCalf must be positive numbers.")
    r = cThigh / cCalf
    return round(r, dp)


def AT_FWR(cForearm: float, cWrist: float, dp: int = 4) -> float:
    """
    Calculate Forearm-to-Wrist Ratio.

    Parameters:
        cForearm (float): Forearm circumference in cm.
        cWrist (float): Wrist circumference in cm.
        dp (int): Number of decimal points. Default = 4.

    Returns:
        float: Forearm-to-Wrist Ratio.
    """
    if cForearm <= 0 or cWrist <= 0:
        raise ValueError("cThigh and cCalf must be positive numbers.")
    r = cForearm / cWrist
    return round(r, dp)


def AT_CHtR(cCalf: float, mHeight:float, dp: int = 4) -> float:
    """
    Calculate Calf-to-Height Ratio.

    Parameters:
        cCalf (float): Calf circumference in cm.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 4.

    Returns:
        float: Calf-to-Height Ratio.
    """
    if mHeight <= 0 or cCalf <= 0:
        raise ValueError("cCalf and mHeight must be positive numbers.")
    r = cCalf / mHeight
    return round(r, dp)


def AT_NHtR(cNeck: float, mHeight:float, dp: int = 4) -> float:
    """
    Calculate Neck-to-Height Ratio.

    Parameters:
        cNeck (float): Neck circumference in cm.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 4.

    Returns:
        float: Neck-to-Height Ratio.
    """
    if mHeight <= 0 or cNeck <= 0:
        raise ValueError("cNeck and mHeight must be positive numbers.")
    r = cNeck / mHeight
    return round(r, dp)


def AI_WHR(cWaist: float, cHip: float, dp: int = 2) -> float:
    """
    Calculate Waist-to-Hip Ratio (WHR).

    Parameters:
        cWaist (float): Waist circumference in cm.
        cHip (float): Hip circumference in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Waist-to-hip ratio.
    """
    if cWaist <= 0 or cHip <= 0:
        raise ValueError("cWaist and cHip circumferences must be positive numbers.")
    whr = cWaist / cHip
    return _round(whr, dp)


def AI_WHtR(cWaist: float, mHeight: float, dp: int = 2) -> float:
    """
    Calculate Waist-to-Height Ratio (WHtR).

    Parameters:
        cWaist (float): Waist circumference in cm.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Waist-to-height ratio.
    """
    if cWaist <= 0 or mHeight <= 0:
        raise ValueError("cWaist circumference and mHeight must be positive numbers.")
    whtr = cWaist / mHeight
    return _round(whtr, dp)


def AI_WTR(cWaist: float, cThigh: float, dp: int = 2) -> float:
    """
    Calculate Waist-to-Thigh Ratio (WTR).

    Parameters:
        cWaist (float): Waist circumference in cm.
        cThigh (float): Thigh in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Waist-to-Thigh ratio.
    """
    if cWaist <= 0 or cThigh <= 0:
        raise ValueError("cWaist and cThigh circumferences must be positive numbers.")
    wtr = cWaist / cThigh
    return _round(wtr, dp)


def AI_ABSI(cWaist: float, mWeight: float, mHeight: float, dp: int = 6) -> float:
    """
    Calculate A Body Shape Index (ABSI).

    Reference: Krakauer, N. Y., & Krakauer, J. C. (2012). A new body shape index predicts mortality hazard independently of body mass index. PloS One 7(7): e39504.

    Parameters:
        cWaist (float): Waist circumference in cm.
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 6.

    Returns:
        float: A Body Shape Index (ABSI).
    """
    if cWaist <= 0 or mWeight <= 0 or mHeight <= 0:
        raise ValueError("cWaist, mWeight, and mHeight must be positive numbers.")
    mHeight_m = mHeight / 100
    cWaist_m = cWaist / 100
    bmi = AT_BMI(mWeight=mWeight, mHeight=mHeight)
    absi = cWaist_m / ((bmi ** (2/3)) * (mHeight_m ** 0.5))
    return _round(absi, dp)


def AI_ConicityIndex(cWaist: float, mWeight: float, mHeight: float, dp: int = 4) -> float:
    """
    Calculate Conicity Index (CI).

    Reference: Valdez, R. (1991). A simple model‑based index of abdominal adiposity. Journal of Clinical Epidemiology, 44(9), 955–956. 

    Parameters:
        cWaist (float): Waist circumference in cm.
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 4.

    Returns:
        float: Conicity Index.
    """
    if cWaist <= 0 or mWeight <= 0 or mHeight <= 0:
        raise ValueError("cWaist, mWeight, and mHeight must be positive numbers.")
    mHeight_m = mHeight / 100
    cWaist_m = cWaist / 100
    # CI = waist / (0.109 * sqrt(weight / height))
    ci = cWaist_m / (0.109 * ((mWeight / mHeight_m) ** 0.5))
    return _round(ci, dp)


def AI_WCR(cWaist: float, cCalf: float, dp: int = 4) -> float:
    """
    Calculate Waist-to-Calf Ratio.

    Parameters:
        cWaist (float): Waist circumference in cm.
        cCalf (float): Calf circumference in cm.
        dp (int): Number of decimal points. Default = 4.

    Returns:
        float: Waist-to-Calf Ratio.
    """
    if cWaist <= 0 or cCalf <= 0:
        raise ValueError("cWaist and cCalf must be positive numbers.")
    r = cWaist / cCalf
    return round(r, dp)


def AI_WHtPI(cWaist: float, mHeight: float, power: float = 2.0, dp: int = 4) -> float:
    """
    Calculate Waist-to-Height Power Index (WHtPI).

    Parameters:
        cWaist (float): Waist circumference in cm.
        mHeight (float): Height in cm.
        power (float): Exponent applied to height (default = 2.0).
        dp (int): Decimal precision. Default = 4.

    Returns:
        float: Waist-to-Height Power Index.
    """
    if mHeight <= 0:
        raise ValueError("mHeight must be positive.")
    if cWaist <= 0:
        raise ValueError("cWaist circumference must be positive.")
    whtpi = cWaist / (mHeight ** power)
    return _round(whtpi, dp)


def AI_BRI(cWaist: float, mHeight: float, dp: int = 2) -> float:
    """
    Calculate Body Roundness Index (BRI).

    Parameters:
        cWaist (float): Waist circumference in cm.
        mHeight (float): Height in cm.
        dp (int): Decimal precision. Default = 2.

    Returns:
        float: Body Roundness Index.
    """
    if mHeight <= 0 or cWaist <= 0:
        raise ValueError("Height and waist must be positive.")
    WC = cWaist / 100.0
    H = mHeight / 100.0
    bri = 364.2 - 365.5 * math.sqrt(1 - (WC / (H * math.pi))**2)
    return _round(bri, dp)


def BFc_USNavy(mHeight: float, cNeck: float, cWaist: float, cHip: float, aGender: str = "male", converter: str = "Siri", dp: int = 2) -> float:
    """
    Calculate body fat percentage using the US Navy method.
    
    Reference: Shake CL, Schlichting C, Mooney LW, Callahan AB, Cohen ME. Predicting percent body fat from circumference measurements. Mil Med. 1993 Jan;158(1):26-31. PMID: 8437737.

    Parameters:
        mHeight (float): Height in cm.
        cNeck (float): Neck circumference in cm.
        cWaist (float): Waist circumference in cm.
        cHip (float, optional): Hip circumference in cm (required for females). Default = 0.
        aGender (str): "male" or "female".
        converter (str): Density to percent body fat converter - "Siri" or "Brozek". Default = Siri.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    converter = converter.lower()
    if mHeight <= 0 or cNeck <= 0 or cWaist <= 0 or (aGender == "female" and cHip <= 0):
        raise ValueError("All measurements must be positive numbers.")
    if aGender == "male":
        if cWaist - cNeck <= 0:
            raise ValueError("cWaist must be larger than cNeck for valid calculation.")
        density = 1.0324 - 0.19077 * math.log10(cWaist - cNeck) + 0.15456 * math.log10(mHeight)
    elif aGender == "female":
        if cWaist + cHip - cNeck <= 0:
            raise ValueError("cWaist + cHip must be larger than cNeck for valid calculation.")
        density = 1.29579 - 0.35004 * math.log10(cWaist + cHip - cNeck) + 0.22100 * math.log10(mHeight)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    if converter == "siri":
        percent_bf = AT_Siri(density)
    else:
        percent_bf = AT_Brozek(density)
    return _round(percent_bf, dp)


def BFc_YMCA(cWaist: float, mWeight: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate body fat percentage using the YMCA method.
    
    Reference: Golding, L.A., Myers, C.R., & Sinning, W.E. (1989). Y’s Way to Physical Fitness: The Complete Guide to Fitness Testing and Instruction. Human Kinetics.

    Parameters:
        cWaist (float): Waist circumference in cm.
        mWeight (float): Body Weight in kg.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    if cWaist <= 0 or mWeight <= 0:
        raise ValueError("cWaist and mWeight must be positive numbers.")
    # Convert to imperial units
    mWeight_lb = mWeight * 2.20462   # kg → lb
    cWaist_in = cWaist * 0.393701    # cm → in
    if aGender == "male":
        percent_bf = 100 * (((4.15 * cWaist_in) - (0.082 * mWeight_lb) - 98.42) / mWeight_lb)
    elif aGender == "female":
        percent_bf = 100 * (((4.15 * cWaist_in) - (0.082 * mWeight_lb) - 76.76) / mWeight_lb)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(percent_bf, dp)


def BFc_mYMCA(cWaist: float, mWeight: float, cWrist: float, cHip: float, cForearm: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate body fat percentage using the modified YMCA method.
    
    Reference: Golding, L.A., Myers, C.R., & Sinning, W.E. (1989). Y’s Way to Physical Fitness: The Complete Guide to Fitness Testing and Instruction. Human Kinetics.

    Parameters:
        cWaist (float): Waist circumference in cm.
        mWeight (float): Body Weight in kg.
        aGender (str): "male" or "female".
        cWrist (float, optional): Wrist circumference in cm (required for females).
        cHip (float, optional): Hip circumference in cm (required for females).
        cForearm (float, optional): Forearm circumference in cm (required for females).
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    if cWaist <= 0 or mWeight <= 0:
        raise ValueError("cWaist and mWeight must be positive numbers.")
    # Convert to imperial units
    mWeight_lb = mWeight * 2.20462   # kg → lb
    cWaist_in = cWaist * 0.393701    # cm → in
    if aGender == "male":
        percent_bf = 100 * (((4.15 * cWaist_in) - (0.082 * mWeight_lb) - 94.42) / mWeight_lb)
    elif aGender == "female":
        if cWrist <= 0 or cHip <= 0 or cForearm <= 0:
            raise ValueError("cWrist, cHip, and cForearm must be positive numbers for females.")
        cWrist_in = cWrist * 0.393701
        cHip_in = cHip * 0.393701
        cForearm_in = cForearm * 0.393701
        percent_bf = 100 * (((0.268 * mWeight_lb) - (0.318 * cWrist_in) + (0.157 * cWaist_in) +
                             (0.245 * cHip_in) - (0.434 * cForearm_in) - 8.987) / mWeight_lb)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(percent_bf, dp)


def BFc_CovertBailey(aAge: int, cWaist: float, cHip: float, cForearm: float, cWrist: float, cThigh: float, cCalf: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate body fat using the Covert Bailey method.
    
    Reference: Bailey, C. (1978). Fit or Fat. Boston: Houghton Mifflin.

    Parameters:
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        cWaist (float, optional): Waist circumference in cm (required for males).
        cHip (float): Hip circumference in cm.
        cForearm (float, optional): Forearm circumference in cm (required for males).
        cWrist (float): Wrist circumference in cm.
        cThigh (float, optional): Thigh circumference in cm (required for females).
        cCalf (float, optional): Calf circumference in cm (required for females).
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat value (method-specific).
    """
    aGender = aGender.lower()
    if aAge <= 0 or cWrist <= 0 or cHip <= 0:
        raise ValueError("aAge, cWrist, and cHip must be positive numbers.")
    # Convert to inches
    cWrist_in = cWrist * 0.393701
    cHip_in = cHip * 0.393701
    if aGender == "male":
        if cWaist <= 0 or cForearm <= 0:
            raise ValueError("cWaist and cForearm must be positive numbers for males.")
        cWaist_in = cWaist * 0.393701
        cForearm_in = cForearm * 0.393701
        if aAge < 30.01:
            r = cWaist_in + (0.5 * cHip_in) - (3.0 * cForearm_in) - cWrist_in
        else:
            r = cWaist_in + (0.5 * cHip_in) - (2.7 * cForearm_in) - cWrist_in
    elif aGender == "female":
        if cThigh <= 0 or cCalf <= 0:
            raise ValueError("cThigh and cCalf must be positive numbers for females.")
        cThigh_in = cThigh * 0.393701
        cCalf_in = cCalf * 0.393701
        if aAge < 30.01:
            r = cHip_in + (0.8 * cThigh_in) - (2.0 * cCalf_in) - cWrist_in
        else:
            r = cHip_in + cThigh_in - (2.0 * cCalf_in) - cWrist_in
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(r, dp)


def BFc_BehnkeWilmore(mWeight: float, cWaist: float, dp: int = 2) -> float:
    """
    Calculate body fat percentage using the Behnke-Wilmore circumference method.
    
    Reference: Wilmore, J.H., & Behnke, A.R. (1969). An anthropometric estimation of body density and lean body weight in young men. Journal of Applied Physiology, 27(1), 25–31. https://doi.org/10.1152/jappl.1969.27.1.25

    Parameters:
        mWeight (float): Body meight in kg.
        cWaist (float): Waist circumference in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    if mWeight <= 0 or cWaist <= 0:
        raise ValueError("mWeight and cWaist must be positive numbers.")
    # Behnke-Wilmore formula (using metric units directly)
    LBW = (1.0817 * mWeight) - (0.7396 * cWaist) + 44.636
    FBW = mWeight - LBW
    bodyfat = (FBW * 100) / mWeight
    return _round(bodyfat, dp)


def BFc_RFM(mHeight: float, cWaist: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate body fat percentage using the Relative Fat Mass (RFM) equation.

    Reference: Woolcott, O.O., Bergman, R.N. (2018) Relative fat mass (RFM) as a new estimator of whole-body fat percentage ─ A cross-sectional study in American adult individuals. Sci Rep 8: 10980. https://doi.org/10.1038/s41598-018-29362-1

    Parameters:
        mHeight (float): Height in cm.
        cWaist (float): Waist circumference in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    if mHeight <= 0 or cWaist <= 0:
        raise ValueError("mHeight and cWaist must be positive numbers.")
    if aGender == "male":
        bf = 64 - (20 * (mHeight / cWaist))
    elif gender == "female":
        bf = 76 - (20 * (mHeight / cWaist))
    else:
        raise ValueError("Gender must be 'male' or 'female'.")
    return _round(bf, dp)


def BFc_BAI(mHeight: float, cHip: float, dp: int = 2) -> float:
    """
    Estimate Body Adiposity Index (BAI).

    Reference: Bergman RN, Stefanovski D, Buchanan TA, Sumner AE, Reynolds JC, Sebring NG, Xiang AH, Watanabe RM. (2011) A better index of body adiposity. Obesity (Silver Spring) 19(5): 1083-9. doi: 10.1038/oby.2011.38. 

    Parameters:
        mHeight (float): Height in cm.
        cHip(float): Hip circumference in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body adiposity index (% body fat).
    """
    if cHip <= 0 or mHeight <= 0:
        raise ValueError("cHip and mHeight must be positive numbers.")
    bai = (cHip / ((mHeight / 100) ** 1.5)) - 18
    return _round(bai, dp)


def BFs_JacksonPollock3(sChest: float, sAbdominal: float, sThigh: float, aAge: int, aGender: str = "male", converter: str = "Siri", dp: int = 2) -> float:
    """
    Calculate body fat percentage using the Jackson-Pollock 3-site method.

    Reference: Jackson, A.S., & Pollock, M.L. (1978). Generalized equations for predicting body density of men. British Journal of Nutrition, 40(3), 497–504. https://doi.org/10.1079/BJN19780152

    Parameters:
        sChest (float): Chest skinfold in mm (required for males).
        sAbdominal (float): Abdominal skinfold in mm.
        sThigh (float): Thigh skinfold in mm.
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        converter (str): Density to percent body fat converter - "Siri" or "Brozek". Default = Siri.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage (via Siri equation).
    """
    aGender = aGender.lower()
    converter = converter.lower()
    if sChest <= 0 or sAbdominal <= 0 or sThigh <= 0 or aAge <= 0:
        raise ValueError("Skinfolds and age must be positive numbers.")
    # Sum of skinfolds
    fold = sChest + sAbdominal + sThigh
    if aGender == "male":
        density = (1.10938 - (0.0008267 * fold) + (0.0000016 * fold * fold) - (0.0002574 * aAge))
    elif aGender == "female":
        density = (1.0994921 - (0.0009929 * fold) + (0.0000023 * fold * fold) - (0.0001392 * aAge))
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    if converter == "siri":
        percent_bf = AT_Siri(density)
    else:
        percent_bf = AT_Brozek(density)
    return _round(percent_bf, dp)


def BFs_JacksonPollock4(sAbdominal: float, sTricep: float, sThigh: float, sSuprailiac: float, aAge: int, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate body fat percentage using the Jackson-Pollock 4-site method.

    Reference: Jackson, A.S., & Pollock, M.L. (1978). Generalized equations for predicting body density of men. British Journal of Nutrition, 40(3), 497–504. https://doi.org/10.1079/BJN19780152

    Parameters:
        sAbdominal (float): Abdominal skinfold in mm.
        sTricep (float): Tricep skinfold in mm.
        sThigh (float): Thigh skinfold in mm.
        sSuprailiac (float): Suprailiac skinfold in mm.
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    if sAbdominal <= 0 or sTricep <= 0 or sThigh <= 0 or sSuprailiac <= 0 or aAge <= 0:
        raise ValueError("Skinfolds and aAge must be positive numbers.")
    # Sum of skinfolds
    fold = sAbdominal + sTricep + sThigh + sSuprailiac
    if aGender == "male":
        bodyfat = ((0.29288 * fold) - (0.0005 * fold * fold) + (0.15845 * aAge) - 5.76377)
    elif aGender == "female":
        bodyfat = ((0.41563 * fold) - (0.00112 * fold * fold) + (0.03661 * aAge) - 4.03653)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(bodyfat, dp)


def BFs_JacksonPollock7(mWeight: float, sChest: float, sAbdominal: float, sThigh: float, sTricep: float, sSubscapular: float, sSuprailiac: float, sMidaxillary: float, aAge: int, aGender: str = "male", converter: str = "Siri", dp: int = 2) -> float:
    """
    Estimate body fat percentage using the Jackson and Pollock 7-site skinfold method.

    Reference: Jackson, AS, Pollock, ML. (1978). Generalized equations for predicting body density of men. British Journal of Nutrition 40(3), 497–504.

    Parameters:
        sChest (float): Chest skinfold (mm).
        sAbdominal (float): Abdominal skinfold (mm).
        sThigh (float): Thigh skinfold (mm).
        sTricep (float): Triceps skinfold (mm).
        sSubscapular (float): Subscapular skinfold (mm).
        sSuprailiac (float): Suprailiac skinfold (mm).
        sMidaxillary (float): Midaxillary skinfold (mm).
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        converter (str): Density to percent body fat converter - "Siri" or "Brozek". Default = Siri.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    converter = converter.lower()
    # Validate inputs
    skinfolds = [sChest, sAbdominal, sThigh, sTricep, sSubscapular, sSuprailiac, sMidaxillary]
    if any(sf <= 0 for sf in skinfolds):
        raise ValueError("All skinfold values must be positive numbers.")
    if aAge <= 0:
        raise ValueError("aAge must be a positive number.")
    sum7 = sum(skinfolds)
    # Body density equations
    if aGender == "male":
        density = 1.112 - (0.00043499 * sum7) + (0.00000055 * (sum7 ** 2)) - (0.00028826 * aAge)
    elif aGender == "female":
        density = 1.097 - (0.00046971 * sum7) + (0.00000056 * (sum7 ** 2)) - (0.00012828 * aAge)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    if converter == "siri":
        percent_bf = AT_Siri(density)
    else:
        percent_bf = AT_Brozek(density)
    return _round(percent_bf , dp)


def BFs_BehnkeWilmore(mWeight: float, sAbdominal: float, dp: int = 2) -> float:
    """
    Calculate body fat percentage using the Behnke-Wilmore skinfold method.

    Reference: Wilmore, J.H., & Behnke, A.R. (1969). An anthropometric estimation of body density and lean body weight in young men. Journal of Applied Physiology, 27(1), 25–31. https://doi.org/10.1152/jappl.1969.27.1.25 

    Parameters:
        mWeight (float): Body Weight in kg.
        sAbdominal (float): Abdominal circumference in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    if mWeight <= 0 or sAbdominal <= 0:
        raise ValueError("mWeight and sAbdominal circumference must be positive numbers.")
    # Behnke-Wilmore (S) formula
    LBW = (0.7927 * mWeight) - (0.3676 * sAbdominal) + 10.260
    FBW = mWeight - LBW
    bodyfat = (FBW * 100) / mWeight
    return _round(bodyfat, dp)


def BFs_DurninWomersley(sBicep: float, sTricep: float, sSubscapular: float, sSuprailiac: float, aAge: int, aGender: str = "male", converter: str = "Siri", dp: int = 2) -> float:
    """
    Estimate body fat percentage using the Durnin–Womersley 4-site skinfold method.

    Reference: Durnin, J.V.G.A. and Womersley, J. (1974). Body fat assessed from the total body density and its estimation from skinfold thickness: measurements on 481 men and women aged from 16 to 72 years. British Journal of Nutrition, 32, 77-97.

    Parameters:
        sBicep (float): Bicep skinfold in mm.
        sTricep (float): Tricep skinfold in mm.
        sSubscapular (float): Subscapular skinfold in mm.
        sSuprailiac (float): Suprailiac skinfold in mm.
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        converter (str): Density to percent body fat converter - "Siri" or "Brozek". Default = Siri.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    converter = converter.lower()
    if sBicep <= 0 or sTricep <= 0 or sSubscapular <= 0 or sSuprailiac <= 0 or aAge <= 0:
        raise ValueError("Skinfolds and age must be positive numbers.")
    # Sum of 4 skinfolds
    fold_sum = sBicep + sTricep + sSubscapular + sSuprailiac
    # Log10 of sum of skinfolds
    log_sum = math.log10(fold_sum)
    # Durnin–Womersley regression equations for body density
    if aGender == "male":
        if 17 <= aAge <= 19: density = 1.1620 - (0.0630 * log_sum)
        elif 20 <= aAge <= 29: density = 1.1631 - (0.0632 * log_sum)
        elif 30 <= aAge <= 39: density = 1.1422 - (0.0544 * log_sum)
        elif 40 <= aAge <= 49: density = 1.1620 - (0.0700 * log_sum)
        elif aAge >= 50: density = 1.1715 - (0.0779 * log_sum)
        else: density = 1.1533 - (0.0643 * log_sum)
    elif aGender == "female":
        if 17 <= aAge <= 19: density = 1.1549 - (0.0678 * log_sum)
        elif 20 <= aAge <= 29: density = 1.1599 - (0.0717 * log_sum)
        elif 30 <= aAge <= 39: density = 1.1423 - (0.0632 * log_sum)
        elif 40 <= aAge <= 49: density = 1.1333 - (0.0612 * log_sum)
        elif age >= 50: density = 1.1339 - (0.0645 * log_sum)
        else: density = 1.1369 - (0.0598 * log_sum)
    else:
        raise ValueError("Gender must be 'male' or 'female'.")
    if converter == "siri":
        percent_bf = AT_Siri(density)
    else:
        percent_bf = AT_Brozek(density)
    return _round(percent_bf, dp)


def BFs_Sloan(sTricep: float, sSuprailiac: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate body fat percentage using the Sloan equation (1967).

    Reference: Sloan AW. (1967) Estimation of body fat in young men. J Appl Physiol. 23(3): 311-5. doi: 10.1152/jappl.1967.23.3.311.

    Parameters:
        sTricep (float): Tricep skinfold thickness in mm.
        sSuprailiac (float): Suprailiac skinfold thickness in mm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    if sTricep <= 0 or sSuprailiac <= 0:
        raise ValueError("Skinfold values must be positive numbers.")
    if aGender == "male":
        # Sloan equation for men
        bf = (0.785 * sTricep) + (0.395 * sSuprailiac) + 5.1
    elif aGender == "female":
        # Sloan equation for women
        bf = (0.610 * sTricep) + (0.417 * sSuprailiac) + 6.1
    else:
        raise ValueError("sGender must be 'male' or 'female'.")
    return _round(bf, dp)


def BFs_Parrillo(mWeight: float, sChest: float, sAbdominal: float, sThigh: float, sTricep: float, sSubscapular: float, sSuprailiac: float, sLowerback: float, sCalf: float, sBicep: float, dp: int = 2) -> float:
    """
    Estimate body fat percentage using the Parrillo 9-site skinfold method.

    Reference: Parrillo, J, Greenwood-Robinson, M. (1993) High-Performance Bodybuilding. Putnam Publishing Group.

    Parameters:
        sChest (float): Chest skinfold (mm).
        sAbdominal (float): Abdominal skinfold (mm).
        sThigh (float): Thigh skinfold (mm).
        sTricep (float): Triceps skinfold (mm).
        sSubscapular (float): Subscapular skinfold (mm).
        sSuprailiac (float): Suprailiac skinfold (mm).
        sLowerback (float): Lower back skinfold (mm).
        sCalf (float): Calf skinfold (mm).
        sBicep (float): Biceps skinfold (mm).
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    # Validate inputs
    skinfolds = [sChest, sAbdominal, sThigh, sTricep, sSubscapular, sSuprailiac, sLowerback, sCalf, sBicep]
    if any(sf <= 0 for sf in skinfolds):
        raise ValueError("All skinfold values must be positive numbers.")
    if mWeight <= 0:
        raise ValueError("mWeight must be positive numbers.")
    mWeight_lb = mWeight * 2.20462   # kg → lb
    total_skinfold = sum(skinfolds)
    bf_percent = (total_skinfold * 27) / mWeight_lb
    return _round(bf_percent, dp)


def BFx_BMI(aBMI: float, aAge: int, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate body fat percentage using the BMI-based method, which is also Deurenberg's equation.

    Reference: Deurenberg, P., Weststrate, J.A., & Seidell, J.C. (1991). Body mass index as a measure of body fatness: age‑ and sex‑specific prediction formulas. European Journal of Clinical Nutrition, 45(5), 271–281. 10.1079/bjn19910073.

    Parameters:
        aBMI (float): Body Mass Index.
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage.
    """
    aGender = aGender.lower()
    if aBMI <= 0 or aAge <= 0:
        raise ValueError("BMI and age must be positive numbers.")
    if aGender == "male":
        sex = 1
    elif aGender == "female":
        sex = 0
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    percent_bf = (1.20 * aBMI) + (0.23 * aAge) - (10.8 * sex) - 5.4
    return _round(percent_bf, dp)

    
def BFx_JacksonPollock3Girth(aAge: int, sChest: float, sAbdominal: float, sThigh: float, cWaist: float, cForearm: float, sTricep: float, sSuprailiac: float, cHip: float, aGender: str = "male", converter: str = "Siri", dp: int = 2) -> float:
    """
    Calculate body fat percentage using the Jackson-Pollock 3-site + girth method.

    Reference: Jackson, A.S., & Pollock, M.L. (1978). Generalized equations for predicting body density of men. British Journal of Nutrition, 40(3), 497–504. https://doi.org/10.1079/BJN19780152

    Parameters:
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        sChest (float, optional): Chest skinfold in mm (required for males).
        sAbdominal (float): Abdominal skinfold in mm.
        sThigh (float): Thigh skinfold in mm.
        cWaist (float, optional): Waist circumference in cm (required for males).
        cForearm (float, optional): Forearm circumference in cm (required for males).
        sTricep (float, optional): Tricep skinfold in mm (required for females).
        sSuprailiac (float, optional): Suprailiac skinfold in mm (required for females).
        cHip (float, optional): Hip circumference in cm (required for females).
        converter (str): Density to percent body fat converter - "Siri" or "Brozek". Default = Siri.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage (via Siri equation).
    """
    aGender = aGender.lower()
    converter = converter.lower()
    if aAge <= 0:
        raise ValueError("Age must be a positive number.")
    if aGender == "male":
        if sChest <= 0 or sAbdominal <= 0 or sThigh <= 0 or cWaist <= 0 or cForearm <= 0:
            raise ValueError("sChest, sAbdominal, cThigh, cWaist, and cForearm must be positive numbers for males.")
        fold = sChest + sAbdominal + sThigh
        density = (1.099075 - (0.0008290 * fold) + (0.0000026 * fold * fold) - (0.0002017 * aAge) - (0.00005675 * cWaist) + (0.00018586 * cForearm))
    elif aGender == "female":
        if sTricep <= 0 or sSuprailiac <= 0 or sThigh <= 0 or cHip <= 0:
            raise ValueError("sTricep, sSuprailiac, cThigh, and cHip must be positive numbers for females.")
        fold = sTricep + sSuprailiac + sThigh
        density = ( 1.1470292 - (0.0009376 * fold) + (0.000003 * fold * fold) - (0.0001392 * aAge) + (0.0005839 * cHip))
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    if converter == "siri":
        percent_bf = AT_Siri(density)
    else:
        percent_bf = AT_Brozek(density)
    return _round(percent_bf, dp)


def BFx_Henry2018(aAge: int, mHeight: float, sBicep: float, sTricep: float, cWaist: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate body fat percentage using the Henry et al. (2018) method for Asian Chinese.

    Reference: Henry, C. J., Shalini, D., Ponnalagu, O., Bi, X., & Tan, S. Y. (2018). New equations to predict body fat in Asian-Chinese adults using age, height, skinfold thickness, and waist circumference. Journal of the Academy of Nutrition and Dietetics, 118(7), 1263-1269.

    Parameters:
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        mHeight (float): Height in cm.
        sBicep (float): Bicep skinfold in mm.
        sTricep (float): Tricep skinfold in mm.
        cWaist (float): Waist circumference in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body fat percentage (via Siri equation).
    """
    aGender = aGender.lower()
    if mHeight <= 0 or sBicep <= 0 or sTricep <= 0 or cWaist <= 0:
            raise ValueError("All measurements must be positive numbers.")
    if aAge <= 0:
        raise ValueError("Age must be a positive number.")
    if aGender == "male":
        percent_bf = (17.81 * math.log10(sTricep + sBicep)) + (0.25 * cWaist) - (0.15 * mHeight) + (0.029 * aAge)
    elif aGender == "female":
        percent_bf = (21.75 * math.log10(sTricep + sBicep)) + (0.22 * cWaist) - (0.099 * mHeight) + (0.064 * aAge)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(percent_bf, dp)


def BMR_MifflinStJeor(mWeight: float, mHeight: float, aAge: int, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate basal metabolic rate (BMR) using the Mifflin-St Jeor equation.

    Reference: Mifflin, M.D., St Jeor, S.T., Hill, L.A., Scott, B.J., Daugherty, S.A., & Koh, Y.O. (1990). A new predictive equation for resting energy expenditure in healthy individuals. American Journal of Clinical Nutrition, 51(2), 241–247. https://doi.org/10.1093/ajcn/51.2.241

    Parameters:
        mWeight (float): Body Weight in kg.
        mHeight (float): Height in cm.
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated BMR (kcal/day).
    """
    aGender = aGender.lower()
    if mWeight <= 0 or mHeight <= 0 or aAge <= 0:
        raise ValueError("mWeight, mHeight, and aAge must be positive numbers.")
    if aGender == "male":
        r = (10 * mWeight) + (6.25 * mHeight) - (5 * aAge) + 5
    elif aGender == "female":
        r = (10 * mWeight) + (6.25 * mHeight) - (5 * aAge) - 161
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(r, dp)


def BMR_HarrisBenedict(mWeight: float, mHeight: float, aAge: int, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate basal metabolic rate (BMR) using the Harris-Benedict equation.

    Reference: Roza, A.M., & Shizgal, H.M. (1984). The Harris Benedict equation reevaluated: resting energy requirements and the body cell mass. American Journal of Clinical Nutrition, 40(1), 168–182. https://doi.org/10.1093/ajcn/40.1.168

    Parameters:
        mWeight (float): Body Weight in kg.
        mHeight (float): Height in cm.
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated BMR (kcal/day).
    """
    aGender = aGender.lower()
    if mWeight <= 0 or mHeight <= 0 or aAge <= 0:
        raise ValueError("mWeight, mHeight, and aAge must be positive numbers.")
    if aGender == "male":
        r = 66.5 + (13.75 * mWeight) + (5.003 * mHeight) - (6.755 * aAge)
    elif aGender == "female":
        r = 655.1 + (9.563 * mWeight) + (1.850 * mHeight) - (4.676 * aAge)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(r, dp)


def BMR_Kleiber(mWeight: float, dp: int = 2) -> float:
    """
    Calculate basal metabolic rate (BMR) using Kleiber's law.

    Reference: Kleiber, M. (1932) Body Size and Metabolism. Hilgardia, 6, 315-353. https://doi.org/10.3733/hilg.v06n11p315

    Parameters:
        mWeight (float): Body Weight in kg.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated BMR (kcal/day).
    """
    if mWeight <= 0:
        raise ValueError("mWeight must be a positive number.")
    # Kleiber's law: BMR ≈ 70 * (mWeight^0.75)
    r = 70 * (mWeight ** 0.75)
    return _round(r, dp)


def BMR_RozaShizgal(mWeight: float, mHeight: float, aAge: int, aGender: str = "male", dp: int = 2) -> float:
    """
    Calculate basal metabolic rate (BMR) using the Roza-Shizgal revision of the Harris-Benedict equation.

    Reference: Roza, A.M., & Shizgal, H.M. (1984). The Harris Benedict equation reevaluated: resting energy requirements and the body cell mass. American Journal of Clinical Nutrition, 40(1), 168–182. https://doi.org/10.1093/ajcn/40.1.168

    Parameters:
        mWeight (float): Body Weight in kg.
        mHeight (float): Height in cm.
        aAge (int): Age in years.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated BMR (kcal/day).
    """
    aGender = aGender.lower()
    if mWeight <= 0 or mHeight <= 0 or aAge <= 0:
        raise ValueError("mWeight, mHeight, and aAge must be positive numbers.")
    if aGender == "male":
        r = 88.362 + (13.397 * mWeight) + (4.799 * mHeight) - (5.677 * aAge)
    elif aGender == "female":
        r = 447.593 + (9.247 * mWeight) + (3.098 * mHeight) - (4.330 * aAge)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(r, dp)


def BMR_Schofield(mWeight: float, aAge: int, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate Basal Metabolic Rate (BMR) using the Schofield (1985) equations.

    Reference: Schofield WN. (1985) Predicting basal metabolic rate, new standards and review of previous work. Hum Nutr Clin Nutr. 39 (Suppl 1): 5-41. 

    Parameters:
        mWeight (float): Body weight in kg.
        aAge (int): age in years.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated BMR (kcal/day).
    """
    aGender = aGender.lower()
    if mWeight <= 0 or aAge <= 0:
        raise ValueError("mWeight and aAge must be positive numbers.")
    # Schofield equations by aAge and sex (mWeight-based)
    if aGender == "male":
        if 3 <= aAge <= 10:
            bmr = (22.706 * mWeight) + 504.3
        elif 11 <= aAge <= 17:
            bmr = (17.686 * mWeight) + 658.2
        elif 18 <= aAge <= 29:
            bmr = (15.057 * mWeight) + 692.2
        elif 30 <= aAge <= 59:
            bmr = (11.472 * mWeight) + 873.1
        elif aAge >= 60:
            bmr = (11.711 * mWeight) + 587.7
        else:
            bmr = (59.512 * mWeight) - 30.4
    elif aGender == "female":
        if 3 <= aAge <= 10:
            bmr = (20.315 * mWeight) + 485.9
        elif 11 <= aAge <= 17:
            bmr = (13.384 * mWeight) + 692.6
        elif 18 <= aAge <= 29:
            bmr = (14.818 * mWeight) + 486.6
        elif 30 <= aAge <= 59:
            bmr = (8.126 * mWeight) + 845.6
        elif aAge >= 60:
            bmr = (9.082 * mWeight) + 658.5
        else:
            bmr = (58.317 * mWeight) - 31.1
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(bmr, dp)


def BMR_Cunningham(mWeight: float, mHeight: float, aGender: str = "male", converter: str = "Boer", dp: int = 2) -> float:
    """
    Estimate Basal Metabolic Rate (BMR) using the Cunningham (1980) equation.

    Reference: Cunningham JJ. (1980) A reanalysis of the factors influencing basal metabolic rate in normal adults. Am J Clin Nutr. 33(11): 2372-4. doi: 10.1093/ajcn/33.11.2372. 

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        converter (str): Lean Body Mass converter - "Boer", "James", or "Hume". Default = Boer.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated BMR (kcal/day).
    """
    aGender = aGender.lower()
    converter = converter.lower()
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mheight must be positive numbers.")
    if converter == "boer":
        lbm = LBM_Boer(mWeight=mWeight, mHeight=mHeight, aGender=aGender)
    elif converter == "james":
        lbm = LBM_James(mWeight=mWeight, mHeight=mHeight, aGender=aGender)
    elif converter == "hume":
        lbm = LBM_Hume(mWeight=mWeight, mHeight=mHeight, aGender=aGender)
    # Cunningham equation
    bmr = 500 + (22 * lbm)
    return _round(bmr, dp)


def BMR_KatchMcArdle(mWeight: float, mHeight: float, aGender: str = "male", converter: str = "Boer", dp: int = 2) -> float:
    """
    Estimate Basal Metabolic Rate (BMR) using the Katch and McArdle (1983) equation.

    Reference: Katch, F., McArdle, W. (1983) Nutrition, Weight Control, and Exercise, 2d. Lea & Febiger.

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        converter (str): Lean Body Mass converter - "Boer", "James", or "Hume". Default = Boer.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated BMR (kcal/day).
    """
    aGender = aGender.lower()
    converter = converter.lower()
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mheight must be positive numbers.")
    if converter == "boer":
        lbm = LBM_Boer(mWeight=mWeight, mHeight=mHeight, aGender=aGender)
    elif converter == "james":
        lbm = LBM_James(mWeight=mWeight, mHeight=mHeight, aGender=aGender)
    elif converter == "hume":
        lbm = LBM_Hume(mWeight=mWeight, mHeight=mHeight, aGender=aGender)
    # Katch and McArdle equation
    bmr = 370 + (21.6 * lbm)
    return _round(bmr, dp)


def BMR_FAO2004(mWeight: float, aAge: int, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate Basal Metabolic Rate (BMR) using the FAO/WHO/UNU (2004) equations.

    Reference: FAO/WHO/UNU Expert Consultation (2004). Human energy requirements: Report of a Joint FAO/WHO/UNU Expert Consultation. Food and Nutrition Technical Report Series, FAO, Rome.

    Parameters:
        mWeight (float): Body mWeight in kg.
        aAge (int): aAge in years.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated BMR (kcal/day).
    """
    aGender = aGender.lower()
    if mWeight <= 0 or aAge <= 0:
        raise ValueError("mWeight and aAge must be positive numbers.")
    # FAO/WHO/UNU (2004) equations by aAge and sex
    if aGender == "male":
        if 3 <= aAge <= 10: bmr = 22.7 * mWeight + 495
        elif 10 < aAge <= 18: bmr = 17.5 * mWeight + 651
        elif 18 < aAge <= 30: bmr = 15.3 * mWeight + 679
        elif 30 < aAge <= 60: bmr = 11.6 * mWeight + 879
        elif aAge > 60: bmr = 13.5 * mWeight + 487
        else: raise ValueError("aAge must be at least 3 years for this equation.")
    elif aGender == "female":
        if 3 <= aAge <= 10: bmr = 22.5 * mWeight + 499
        elif 10 < aAge <= 18: bmr = 12.2 * mWeight + 746
        elif 18 < aAge <= 30: bmr = 14.7 * mWeight + 496
        elif 30 < aAge <= 60: bmr = 8.7 * mWeight + 829
        elif aAge > 60: bmr = 10.5 * mWeight + 596
        else: raise ValueError("aAge must be at least 3 years for this equation.")
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(bmr, dp)


def BSA_DuBois(mWeight: float, mHeight: float, dp: int = 2) -> float:
    """
    Estimate body surface area (BSA) using the Du Bois & Du Bois (1916) formula.

    Reference: Du Bois D, Du Bois EF. A formula to estimate the approximate surface area if height and weight be known. 1916. Nutrition. 1989 Sep-Oct;5(5):303-11; discussion 312-3. PMID: 2520314.

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body surface area (m²).
    """
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mheight must be positive numbers.")
    bsa = 0.007184 * (mHeight ** 0.725) * (mWeight ** 0.425)
    return _round(bsa, dp)


def BSA_Mosteller(mWeight: float, mHeight: float, dp: int = 2) -> float:
    """
    Estimate body surface area (BSA) using the Mosteller (1987) formula.

    Reference: Mosteller RD. (1987) Simplified calculation of body-surface area. N Engl J Med 317(17):1098. doi: 10.1056/NEJM198710223171717.

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body surface area (m²).
    """
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mheight must be positive numbers.")
    bsa = ((mHeight * mWeight) / 3600) ** 0.5
    return _round(bsa, dp)


def BSA_Haycock(mWeight: float, mHeight: float, dp: int = 2) -> float:
    """
    Estimate body surface area (BSA) using the Haycock et al. (1978) formula.

    Reference: Haycock GB, Schwartz GJ, Wisotsky DH. (1978) Geometric method for measuring body surface area: a height-weight formula validated in infants, children, and adults. J Pediatr 93(1):62-6. doi: 10.1016/s0022-3476(78)80601-5. 

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated body surface area (m²).
    """
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mheight must be positive numbers.")
    bsa = 0.024265 * (mWeight ** 0.5378) * (mHeight ** 0.3964)
    return _round(bsa, dp)


def BSA_GehanGeorge(mWeight: float, mHeight: float, dp: int = 4) -> float:
    """
    Estimate Body Surface Area (BSA) using the Gehan & George (1970) formula.

    Reference: Gehan, EA, George, SL. (1970) Estimation of human body surface area from height and weight. Cancer Chemotherapy Reports, 54(4), 225–235.

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 4.

    Returns:
        float: Estimated body surface area (m^2).
    """
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mHeight must be positive numbers.")
    # Gehan & George formula
    bsa = 0.0235 * (mHeight ** 0.42246) * (mWeight ** 0.51456)
    return _round(bsa, dp)


def BSA_Fujimoto(mWeight: float, mHeight: float, dp: int = 4) -> float:
    """
    Estimate Body Surface Area (BSA) using the Fujimoto formula.

    Reference: Kouno T, Katsumata N, Mukai H, Ando M, Watanabe T. (2003) Standardization of the body surface area (BSA) formula to calculate the dose of anticancer agents in Japan. Jpn J Clin Oncol.33(6):309-13. doi: 10.1093/jjco/hyg062

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 4.

    Returns:
        float: Estimated body surface area (m^2).
    """
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mHeight must be positive numbers.")
    # Fujimoto formula
    bsa = 0.008883 * (mHeight ** 0.663) * (mWeight ** 0.444)
    return _round(bsa, dp)


def LBM_Boer(mWeight: float, mHeight: float, aGender: str = "male",dp: int = 2) -> float:
    """
    Estimate Lean Body Mass (LBM) using the Boer (1984) formula.

    Reference: Boer P. (1987). Estimated lean body mass as an index for normalization of body fluid volumes in humans. Am J Physiol. 247(4 Pt 2):F632-6. doi: 10.1152/ajprenal.1984.247.4.F632.

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated lean body mass (kg).
    """
    aGender = aGender.lower()
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mheight must be positive numbers.")
    if aGender == "male":
        lbm = (0.407 * mWeight) + (0.267 * mHeight) - 19.2
    elif aGender == "female":
        lbm = (0.252 * mWeight) + (0.473 * mHeight) - 48.3
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(lbm, dp)


def LBM_James(mWeight: float, mHeight: float, aGender: str = "male",dp: int = 2) -> float:
    """
    Estimate Lean Body Mass (LBM) using the James (1977) formula.

    Reference: James, W.P.T. (1977), Research on obesity. Nutrition Bulletin, 4: 187-190. https://doi.org/10.1111/j.1467-3010.1977.tb00966.x

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated lean body mass (kg).
    """
    aGender = aGender.lower()
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mheight must be positive numbers.")

    if aGender == "male":
        lbm = (1.10 * mWeight) - (128 * (mWeight ** 2 / (mHeight ** 2)))
    elif aGender == "female":
        lbm = (1.07 * mWeight) - (148 * (mWeight ** 2 / (mHeight ** 2)))
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(lbm, dp)


def LBM_Hume(mWeight: float, mHeight: float, aGender: str = "male",dp: int = 2) -> float:
    """
    Estimate Lean Body Mass (LBM) using the Hume (1966) formula.

    Reference: Hume, R. (1966). Prediction of lean body mass from height and weight. Journal of Clinical Pathology, 19(4), 389–391. https://doi.org/10.1136/jcp.19.4.389

    Parameters:
        mWeight (float): Body weight in kg.
        mHeight (float): Height in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated lean body mass (kg).
    """
    aGender = aGender.lower()
    if mWeight <= 0 or mHeight <= 0:
        raise ValueError("mWeight and mheight must be positive numbers.")
    if aGender == "male":
        lbm = (0.32810 * mWeight) + (0.33929 * mHeight) - 29.5336
    elif aGender == "female":
        lbm = (0.29569 * mWeight) + (0.41813 * mHeight) - 43.2933
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(lbm, dp)


def IBW_Robinson(mHeight: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate Ideal Body Weight (IBW) using the Robinson (1983) formula.

    Reference: Robinson, J. D., Lupkiewicz, S. M., Palenik, L., Lopez, L. M., & Ariet, M. (1983). Determination of ideal body weight for drug dosage calculations. American Journal of Health-System Pharmacy, 40(6), 1016–1019. https://doi.org/10.1093/AJHP/40.6.1016

    Parameters:
        mHeight (float): Height in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated ideal body weight (kg).
    """
    aGender = aGender.lower()
    if mHeight <= 0:
        raise ValueError("mHeight must be a positive number.")
    # Convert height to inches
    height_in = mHeight * 0.393701
    # Height above 5 feet (60 inches)
    extra_in = height_in - 60
    if aGender == "male":
        iw = 52.0 + (1.9 * extra_in)
    elif aGender == "female":
        iw = 49.0 + (1.7 * extra_in)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(iw, dp)


def IBW_Miller(mHeight: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate Ideal Body Weight (IBW) using the Miller formula.

    Reference: Miller, D. R., Carlson, J. D., Lloyd, B. J., & Day, B. J. (1983). Determining ideal body-weight (and mass). American Journal of Hospital Pharmacy, 40(10), 1622-1622.

    Parameters:
        mHeight (float): Height in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated ideal body weight (kg).
    """
    aGender = aGender.lower()
    if mHeight <= 0:
        raise ValueError("mHeight must be a positive number.")
    # Convert height to inches
    height_in = mHeight * 0.393701
    # Height above 5 feet (60 inches)
    extra_in = height_in - 60
    if aGender == "male":
        iw = 56.2 + (1.41 * extra_in)
    elif aGender == "female":
        iw = 53.1 + (1.36 * extra_in)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(iw, dp)


def IBW_Hamwi(mHeight: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate Ideal Body Weight (IBW) using the Hamwi formula.

    Reference: Hamwi, G. J. (1964). Changing dietary concepts. Diabetes mellitus: diagnosis and treatment, 1, 73-8.

    Parameters:
        mHeight (float): Height in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated ideal body weight (kg).
    """
    aGender = aGender.lower()
    if mHeight <= 0:
        raise ValueError("mHeight must be a positive number.")
    # Convert height to inches
    height_in = mHeight * 0.393701
    # Height above 5 feet (60 inches)
    extra_in = height_in - 60
    if aGender == "male":
        iw = 48.0 + (2.7 * extra_in)
    elif aGender == "female":
        iw = 45.5 + (2.2 * extra_in)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(iw, dp)


def IBW_Devine(mHeight: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate Ideal Body Weight (IBW) using the Devine formula.

    Reference: Devine, BJ. (1974) Gentamicin therapy. Drug Intelligence & Clinical Pharmacy 8: 650–655.

    Parameters:
        mHeight (float): Height in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated ideal body weight (kg).
    """
    aGender = aGender.lower()
    if mHeight <= 0:
        raise ValueError("mHeight must be a positive number.")
    # Convert height to inches
    height_in = mHeight * 0.393701
    # Height above 5 feet (60 inches)
    extra_in = height_in - 60
    if aGender == "male":
        iw = 50.0 + (2.3 * extra_in)
    elif aGender == "female":
        iw = 45.5 + (2.3 * extra_in)
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(iw, dp)


def IBW_Lemmens(mHeight: float, dp: int = 2) -> float:
    """
    Estimate Ideal Body Weight (IBW) using the Lemmens et al. (2005) formula.

    Reference: Lemmens HJ, Brodsky JB, Bernstein DP. (2005) Estimating ideal body weight--a new formula. Obes Surg. 15(7): 1082-3. doi: 10.1381/0960892054621350. 

    Parameters:
        mHeight (float): Height in cm.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated ideal body weight (kg).
    """
    if mHeight <= 0:
        raise ValueError("mHeight must be a positive number.")
    iw = 22 * (mHeight / 100) * (mHeight / 100)
    return _round(iw, dp)


def IBW_Peterson(mHeight: float, mWeight: float, dp: int = 2) -> float:
    """
    Estimate Ideal Body Weight (IBW) using the Peterson et al. (2016) formula.

    Reference: Peterson, C. M., Thomas, D. M., Blackburn, G. L., & Heymsfield, S. B. (2016). Universal equation for estimating ideal body weight and body weight at any BMI. The American journal of clinical nutrition, 103(5), 1197-1203.

    Parameters:
        mHeight (float): Height in cm.
        mWeight (float): Body weight in kg.
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated ideal body weight (kg).
    """
    if mHeight <= 0 or mWeight <= 0:
        raise ValueError("mHeight and mWeight must be positive numbers.")
    BMI = AT_BMI(mHeight=mHeight, mWeight=mWeight)
    iw = (2.2 * BMI) + (3.5 * BMI * ((mHeight / 100) - 1.5))
    return _round(iw, dp)


def IBW_KeysBrozek(mHeight: float, aGender: str = "male", dp: int = 2) -> float:
    """
    Estimate Ideal Body Weight (IBW) using the Keys & Brozek (1953) formula.

    Reference: Keys, A, Brozek, J. (1953). Body fat in adult man. Physiological Reviews, 33(3), 245–325. https://doi.org/10.1152/physrev.1953.33.3.245

    Parameters:
        mHeight (float): Height in cm.
        aGender (str): "male" or "female".
        dp (int): Number of decimal points. Default = 2.

    Returns:
        float: Estimated ideal body weight (kg).
    """
    aGender = aGender.lower()
    if mHeight <= 0:
        raise ValueError("mHeight must be a positive number.")
    if aGender == "male":
        iw = (0.73 * mHeight) - 59.4
    elif aGender == "female":
        iw = (0.65 * mHeight) - 50.0
    else:
        raise ValueError("aGender must be 'male' or 'female'.")
    return _round(iw, dp)


def SSF(sAbdominal: float, sBicep: float, sCalf: float, sChest: float, sLowerback: float, sMidaxillary: float, sSubscapular: float, sSuprailiac: float, sThigh: float, sTricep: float, foldset: str | list = "Sum3_men", dp: int = 2) -> float:
    """
    Calculate Sum of Skinfolds (SSF) based on predefined or custom fold sets.

    Parameters:
        sAbdominal, sBicep, sCalf, sChest, sLowerback, sMidaxillary,
        sSubscapular, sSuprailiac, sThigh, sTricep (float): Skinfolds in mm.
        foldset (str | list): Either a predefined set name ("Sum2", "Sum3_men", "Sum3_women", "Sum4", "Sum5", "Sum7", "Sum9", "All") or a custom list of site names.
        dp (int): Decimal precision. Default = 2.

    Returns:
        float: Sum of selected skinfolds (mm).
    """
    sets = {
        "Sum2": ["sTricep", "sSubscapular"],
        "Sum3_men": ["sChest", "sAbdominal", "sThigh"],  # Jackson–Pollock men
        "Sum3_women": ["sTricep", "sSuprailiac", "sThigh"],  # Jackson–Pollock women
        "Sum4": ["sTricep", "sBicep", "sSubscapular", "sSuprailiac"],  # Durnin–Womersley
        "Sum5": ["sTricep", "sSubscapular", "sSuprailiac", "sAbdominal", "sThigh"],
        "Sum7": ["sChest", "sMidaxillary", "sTricep", "sSubscapular",
                 "sAbdominal", "sSuprailiac", "sThigh"],  # JP7
        "Sum9": ["sChest", "sAbdominal", "sThigh", "sTricep", "sSubscapular",
                 "sSuprailiac", "sLowerback", "sCalf", "sBicep"],  # Parrillo
        "All": ["sAbdominal", "sBicep", "sCalf", "sChest", "sLowerback",
                "sMidaxillary", "sSubscapular", "sSuprailiac", "sThigh", "sTricep"]
    }
    # Map argument names to values
    values = {
        "sAbdominal": sAbdominal,
        "sBicep": sBicep,
        "sCalf": sCalf,
        "sChest": sChest,
        "sLowerback": sLowerback,
        "sMidaxillary": sMidaxillary,
        "sSubscapular": sSubscapular,
        "sSuprailiac": sSuprailiac,
        "sThigh": sThigh,
        "sTricep": sTricep
    }
    # Resolve foldset
    if isinstance(foldset, str):
        if foldset not in sets:
            raise ValueError(f"Unknown foldset '{foldset}'. Choose from {list(sets.keys())}")
        sites = sets[foldset]
    elif isinstance(foldset, list):
        sites = foldset
        for site in sites:
            if site not in values:
                raise ValueError(f"Unknown site '{site}' in custom foldset.")
    else:
        raise TypeError("foldset must be a string or a list of site names.")
    total = sum(values[site] for site in sites)
    return round(total, dp)


def FC_Metropolitan(mHeight: float, cWrist: float, aGender: str) -> str:
    """
    Classify body frame size using Body Frame Index (Metropolitan Life, 1983).
    
    Reference: Insurance, M. L. (1983). Metropolitan height and weight tables. Statistical Bulletin, 64, 2-9.

    Parameters:
        mHeight (float): Height in centimeters.
        cWrist (float): Wrist circumference in centimeters.
        aGender (str): "male" or "female".

    Returns:
        str: "Small", "Medium", or "Large" frame size.
    """
    aGender = aGender.lower()
    if mHeight <= 0 or cWrist <= 0:
        raise ValueError("Height and wrist circumference must be positive values.")
    r = mHeight / cWrist
    if aGender == "male":
        if r > 10.4: return "Small"
        elif 9.6 <= r <= 10.4: return "Medium"
        else: return "Large"
    elif aGender == "female":
        if r > 11.0: return "Small"
        elif 10.1 <= r <= 11.0: return "Medium"
        else: return "Large"
    else:
        return "Undefined"


def ST_HeathCarter(sTricep: float, sSubscapular: float, sSuprailiac: float, bHumerus: float, bFemur: float,cBicep: float, cCalf: float, mHeight: float, mWeight: float, dp: int = 2) -> dict:
    """
    Classify somatotype using Heath–Carter method.

    Reference: Carter, J. L., & Heath, B. H. (1990). Somatotyping: development and applications (Vol. 5). Cambridge university press.

    Parameters:
        sTricep (float): Triceps skinfold (mm).
        sSubscapular (float): Subscapular skinfold (mm).
        sSuprailiac (float): Suprailiac skinfold (mm).
        bHumerus (float): Humerus bicondylar breadth (cm).
        bFemur (float): Femur bicondylar breadth (cm).
        cBicep (float): Flexed arm girth (cm).
        cCalf (float): Calf girth (cm).
        mHeight (float): Height (cm).
        mWeight (float): Body weight (kg).
        dp (int): Decimal precision. Default = 2.

    Returns:
        dict: Endomorphy, Mesomorphy, Ectomorphy values.
    """
    # Endomorphy (fatness)
    X = sTricep + sSubscapular + sSuprailiac
    endo = -0.7182 + (0.1451 * X) - (0.00068 * (X**2)) + (0.0000014 * (X**3))
    # Mesomorphy (musculoskeletal robustness)
    meso = (0.858 * bHumerus) + (0.601 * bFemur) + (0.188 * cBicep) + (0.161 * cCalf) - (0.131 * (mHeight/10)) + 4.5
    # Ectomorphy (linearity)
    HWR = mHeight / (mWeight**(1/3))
    if HWR >= 40.75: ecto = 0.732*HWR - 28.58
    elif HWR > 38.25: ecto = 0.463*HWR - 17.63
    else: ecto = 0.1
    return {
        "Endomorphy": _round(endo, dp),
        "Mesomorphy": _round(meso, dp),
        "Ectomorphy": _round(ecto, dp)
    }


measurements = {"aAge": None,"aBMI": None, "aDensity": None, "aGender": None, "mHeight": None, "mWeight": None, "cBicep": None, "cCalf": None, "cForearm": None, "cHip": None, "cNeck": None, "cThigh": None, "cWaist": None, "cWrist": None, "bFemur": None, "bHumerus": None, "sAbdominal": None, "sBicep": None, "sCalf": None, "sChest": None, "sLowerback": None, "sMidaxillary": None, "sSubscapular": None, "sSuprailiac": None, "sThigh": None, "sTricep": None}


def tester():
    measurements["aGender"] = "male"   # gender
    measurements["aAge"] = 36          # years
    measurements["mHeight"] = 173      # cm
    measurements["mWeight"] = 75       # kg
    measurements["cNeck"] = 38         # cm
    measurements["cBicep"] = 29        # cm
    measurements["cForearm"] = 25      # cm
    measurements["cWrist"] = 17        # cm
    measurements["cWaist"] = 85        # cm
    measurements["cHip"] = 95          # cm
    measurements["cThigh"] = 53        # cm
    measurements["cCalf"] = 38         # cm
    measurements["bHumerus"] = 6.8     # cm
    measurements["bFemur"] = 9.0       # cm
    measurements["sChest"] = 12        # mm skinfold
    measurements["sAbdominal"] = 15    # mm skinfold
    measurements["sThigh"] = 12        # mm skinfold
    measurements["sCalf"] = 8          # mm skinfold 
    measurements["sBicep"] = 8         # mm skinfold 
    measurements["sTricep"] = 13       # mm skinfold 
    measurements["sMidaxillary"] = 14  # mm skinfold 
    measurements["sSubscapular"] = 14  # mm skinfold 
    measurements["sSuprailiac"] = 14   # mm skinfold 
    measurements["sLowerback"] = 14    # mm skinfold 
    results = profiler(measurements)
    report(results)


def report(results):
    """
    Nicely print profiler results grouped by category.
    Accepts a dictionary of metrics and prints them in sections.
    """
    # Define categories and group keys
    categories = {
        "Anthropometric Indices": ["AT_BMI", "AT_PI", "AT_TCR", "AT_FWR", "AT_CHtR", "AT_NHtR"],
        "Adiposity Indices": ["AI_WHR", "AI_WHtR", "AI_WTR", "AI_ABSI", "AI_ConicityIndex",
                              "AI_WCR", "AI_WHtPI", "AI_BRI"],
        "Body Fat (Circumference)": ["BFc_USNavy", "BFc_YMCA", "BFc_mYMCA", "BFc_CovertBailey",
                                     "BFc_BehnkeWilmore", "BFc_RFM", "BFc_BAI"],
        "Body Fat (Skinfold)": ["BFs_JacksonPollock3", "BFs_JacksonPollock4", "BFs_JacksonPollock7",
                                "BFs_BehnkeWilmore", "BFs_DurninWomersley", "BFs_Sloan", "BFs_Parrillo"],
        "Body Fat (Other)": ["BFx_BMI", "BFx_JacksonPollock3Girth", "BFx_Henry2018", "BF-Average"],
        "Basal Metabolic Rate": ["BMR_MifflinStJeor", "BMR_HarrisBenedict", "BMR_Kleiber",
                                 "BMR_RozaShizgal", "BMR_Schofield", "BMR_Cunningham",
                                 "BMR_KatchMcArdle", "BMR_FAO2004", "BMR-Average"],
        "Body Surface Area": ["BSA_DuBois", "BSA_Mosteller", "BSA_Haycock",
                              "BSA_GehanGeorge", "BSA_Fujimoto", "BSA-Average"],
        "Lean Body Mass": ["LBM_Boer", "LBM_James", "LBM_Hume", "LBM-Average"],
        "Ideal Body Weight": ["IBW_Robinson", "IBW_Miller", "IBW_Hamwi", "IBW_Devine",
                              "IBW_Lemmens", "IBW_Peterson", "IBW_KeysBrozek", "IBW-Average"],
        "Frame Classification": ["FC_Metropolitan"],
        "Somatotype (Heath-Carter)": ["ST_HeathCarter (Endomorphy)",
                                      "ST_HeathCarter (Mesomorphy)",
                                      "ST_HeathCarter (Ectomorphy)"]
    }
    units = {
        "Anthropometric Indices": {"AT_BMI": "kg/m²", "AT_PI": "kg/m³", "AT_TCR": "", "AT_FWR": "", "AT_CHtR": "", "AT_NHtR": ""},
        "Adiposity Indices": {"AI_WHR": "", "AI_WHtR": "", "AI_WTR": "", "AI_ABSI": "", "AI_ConicityIndex": "", "AI_WCR": "", "AI_WHtPI": "", "AI_BRI": ""},
        "Body Fat (Circumference)": {"BFc_USNavy": "%", "BFc_YMCA": "%", "BFc_mYMCA": "%", "BFc_CovertBailey": "%", "BFc_BehnkeWilmore": "%", "BFc_RFM": "%", "BFc_BAI": "%"},
        "Body Fat (Skinfold)": {"BFs_JacksonPollock3": "%", "BFs_JacksonPollock4": "%", "BFs_JacksonPollock7": "%", "BFs_BehnkeWilmore": "%", "BFs_DurninWomersley": "%", "BFs_Sloan": "%", "BFs_Parrillo": "%"},
        "Body Fat (Other)": {"BFx_BMI": "%", "BFx_JacksonPollock3Girth": "%", "BFx_Henry2018": "%", "BF-Average": "%"},
        "Basal Metabolic Rate": {"BMR_MifflinStJeor": "kcal/day", "BMR_HarrisBenedict": "kcal/day", "BMR_Kleiber": "kcal/day", "BMR_RozaShizgal": "kcal/day", "BMR_Schofield": "kcal/day", "BMR_Cunningham": "kcal / day", "BMR_KatchMcArdle": "kcal/day", "BMR_FAO2004": "kcal/day", "BMR-Average": "kcal/day"},
        "Body Surface Area": {"BSA_DuBois": "m²", "BSA_Mosteller": "m²", "BSA_Haycock": "m²",
                                      "BSA_GehanGeorge": "m²", "BSA_Fujimoto": "m²", "BSA-Average": "m²"},
        "Lean Body Mass": {"LBM_Boer": "kg", "LBM_James": "kg", "LBM_Hume": "kg", "LBM-Average": "kg"},
        "Ideal Body Weight": {"IBW_Robinson": "kg", "IBW_Miller": "kg", "IBW_Hamwi": "kg", "IBW_Devine": "kg", "IBW_Lemmens": "kg", "IBW_Peterson": "kg", "IBW_KeysBrozek": "kg", "IBW-Average": "kg"},
        "Frame Classification": {"FC_Metropolitan": ""},
        "Somatotype (Heath-Carter)": {"ST_HeathCarter (Endomorphy)": "", "ST_HeathCarter (Mesomorphy)": "", "ST_HeathCarter (Ectomorphy)": ""}
        }
    description = {
        "Anthropometric Indices": "(a) AT_BMI: Body Mass Index (BMI). (b) AT_PI: Ponderal Index (PI). (c) AT_TCR: Thigh-to-Calf Ratio. (d) AT_FWR: Forearm-to-Wrist Ratio. (e) AT_CHtR: Calf-to-Height Ratio (f) AT_NHtR: Neck-to-Height Ratio",
        "Adiposity Indices": "(a) AI_WHR: Waist-to-Hip Ratio. (b) AI_WHtR: Waist-to-Height Ratio. (c) AI_WTR: Waist-to-Thigh Ratio. (d) AI_ABSI: A Body Shape Index. (e) AI_ConicityIndex: Conicity Index. (f) AI_WCR: Waist-to-Calf Ratio. (g) AI_WHtPI: Waist-to-Height Power Index. (h) AI_BRI: Body Roundness Index",
        "Body Fat (Circumference)": "(a) BFc_USNavy: US Navy method. (b) BFc_YMCA: YMCA method. (c) BFc_mYMCA: Modified YMCA method. (d) BFc_CovertBailey: Covert Bailey method. (e) BFc_BehnkeWilmore: Behnke-Wilmore circumference method. (f) BFc_RFM: Relative Fat Mass. (g) BFc_BAI: Body Adiposity Index.",
        "Body Fat (Skinfold)": "(a) BFs_JacksonPollock3: Jackson-Pollock 3-site method. (b) BFs_JacksonPollock4: Jackson-Pollock 4-site method. (c) BFs_JacksonPollock7: Jackson-Pollock 7-site method. (d) BFs_BehnkeWilmore: Behnke-Wilmore skinfold method. (e) BFs_DurninWomersley: Durnin–Womersley 4-site skinfold method. (f) BFs_Sloan: Sloan equation (1967). (g) BFs_Parrillo: Parrillo 9-site formula.",
        "Body Fat (Other)": "(a) BFx_BMI: Estimate from BMI using Deurenberg's equation. (b) BFx_JacksonPollock3Girth: Jackson-Pollock 3-site + girth method.  (c) BFx_Henry2018: Henry et al. (2018) method for Asian Chinese.",
        "Basal Metabolic Rate": "(a) BMR_MifflinStJeor: Mifflin-St Jeor equation. (b) BMR_HarrisBenedict: Harris-Benedict equation. (c) BMR_Kleiber: Kleiber's law. (d) BMR_RozaShizgal: Roza-Shizgal revision of the Harris-Benedict equation. (e) BMR_Schofield: Schofield equation. (f) BMR_Cunningham: Cunningham (1980) equation. (g) BMR_KatchMcArdle: Katch and McArdle (1983) equation. (h) BMR_FAO2004: FAO/WHO/UNU (2004) equations.",
        "Body Surface Area": "(a) BSA_DuBois: Du Bois & Du Bois (1916) formula. (b) BSA_Mosteller: Mosteller (1987) formula. (c) BSA_Haycock: Haycock et al. (1978) formula. (d) BSA_GehanGeorge: Gehan & George (1970) formula. (e) BSA_Fujimoto: Fujimoto formula",
         "Lean Body Mass":  "(a) LBM_Boer: Boer (1984) formula. (b) LBM_James: James (1977) formula. (c) LBM_Hume: Hume (1966) formula.",
        "Ideal Body Weight": "(a) IBW_Robinson: Robinson (1983) formula. (b) IBW_Miller: Miller et al. (1983) formula. (c) IBW_Hamwi: Hamwi (1964) formula. (d) IBW_Devine: Devine (1974) formula. (e) IBW_Lemmens: Lemmens et al. (2005) formula. (f) IBW_Peterson: Peterson et al. (2016) formula. (g) IBW_KeysBrozek: Keys & Brozek (1953) formula.",
        "Frame Classification": "(a) FC_Metropolitan: Metropolitan Life (1983)",
        "Somatotype (Heath-Carter)": "(a) ST_HeathCarter: Heath and Carter method."
    }

    # Print report
    print("="*50)
    print("Anthropometrical Report by AnthroPy")
    print("="*50)

    for category, keys in categories.items():
        # Filter only available metrics
        available = {k: results[k] for k in keys if k in results}
        if available:
            print(f"\n{category}:")
            print("-"*len(category))
            for k, v in available.items():
                unit = units[category][k]
                print(f"  {k:35s}: {v} {unit}")
        descriptor = description[category]
        print(f"Description: {descriptor}")
    print("\nEnd of Report")
    print("="*50)


def profiler(measurements):
    """
    Compute anthropometric, body composition, and metabolic metrics based on the provided measurement parameters.

    The function accepts a dictionary of measurement values (e.g., weight, height, waist circumference, age, sex, etc.) and attempts to calculate all possible indices and derived metrics. If the required parameters for a given metric are missing, that metric will remain as None in the results dictionary.

    Parameters
    ----------
    measurements : dict
        Dictionary of measurement inputs. Keys should correspond to expected parameter names (e.g., "mWeight", "mHeight", "cWaist", "cHip", "aGender", "mAge", "cNeck", etc.). Values should be numeric (floats/ints) or categorical (e.g., "M"/"F" for gender). Complete measurement dictionary is {"aAge": None,"aBMI": None, "aDensity": None, "aGender": None, "mHeight": None, "mWeight": None, "cBicep": None, "cCalf": None, "cForearm": None, "cHip": None, "cNeck": None, "cThigh": None, "cWaist": None, "cWrist": None, "bFemur": None, "bHumerus": None, "sAbdominal": None, "sBicep": None, "sCalf": None, "sChest": None, "sLowerback": None, "sMidaxillary": None, "sSubscapular": None, "sSuprailiac": None, "sThigh": None, "sTricep": None}. Missing or None values are ignored.

    Returns
    -------
    results : dict
        Dictionary containing all supported metrics. Keys are metric names (e.g., "AT_BMI", "AI_WHR", "BFc_USNavy", "BMR_MifflinStJeor", "BSA_DuBois", "LBM_Boer", "IBW_Robinson", "ST_HeathCarter (Endomorphy)"). Values are floats, strings, or None if the metric could not be computed due to insufficient input parameters.

    Notes
    -----
    - Metrics are grouped into categories:
        * Anthropometric indices (AT_*)
        * Adiposity indices (AI_*)
        * Body fat estimates (BFc_*, BFs_*, BFx_*)
        * Basal metabolic rate (BMR_*)
        * Body surface area (BSA_*)
        * Lean body mass (LBM_*)
        * Ideal body weight (IBW_*)
        * Frame classification (FC_*)
        * Somatotype (ST_*)
    - Average values (e.g., "BF-Average", "BMR-Average", "BSA-Average", "LBM-Average", "IBW-Average") are computed when multiple methods within a category are available.
    """
    given_parameters = [k for k, v in measurements.items() if v is not None]
    BF, BMR, BSA, LBM, IBW = [], [], [], [], []
    results = {
        "AT_BMI": None, "AT_PI": None, "AT_TCR": None, "AT_FWR": None, "AT_CHtR": None, "AT_NHtR": None, 
        "AI_WHR": None, "AI_WHtR": None, "AI_WTR": None, "AI_ABSI": None, "AI_ConicityIndex": None, "AI_WCR": None, "AI_WHtPI": None, "AI_BRI": None, 
        "BFc_USNavy": None, "BFc_YMCA": None, "BFc_mYMCA": None, "BFc_CovertBailey": None, "BFc_BehnkeWilmore": None, "BFc_RFM": None, "BFc_BAI": None, 
        "BFs_JacksonPollock3": None, "BFs_JacksonPollock4": None, "BFs_JacksonPollock7": None, "BFs_BehnkeWilmore": None, "BFs_DurninWomersley": None, "BFs_Sloan": None, "BFs_Parrillo": None, 
        "BFx_BMI": None, "BFx_JacksonPollock3Girth": None, "BFx_Henry2018": None, 
        "BF-Average": None,
        "BMR_MifflinStJeor": None, "BMR_HarrisBenedict": None, "BMR_Kleiber": None, "BMR_RozaShizgal": None, "BMR_Schofield": None, "BMR_Cunningham": None, "BMR_KatchMcArdle": None, "BMR_FAO2004": None, 
        "BMR-Average": None,
        "BSA_DuBois": None, "BSA_Mosteller": None, "BSA_Haycock": None, "BSA_GehanGeorge": None, "BSA_Fujimoto": None, 
        "BSA-Average": None,
        "LBM_Boer": None, "LBM_James": None, "LBM_Hume": None, 
        "LBM-Average": None,
        "IBW_Robinson": None, "IBW_Miller": None, "IBW_Hamwi": None, "IBW_Devine": None, "IBW_Lemmens": None, "IBW_Peterson": None, "IBW_KeysBrozek": None, 
        "IBW-Average": None,
        "SSF": None, 
        "FC_Metropolitan": None, 
        "ST_HeathCarter (Endomorphy)": None, "ST_HeathCarter (Mesomorphy)": None, "ST_HeathCarter (Ectomorphy)": None, 
        }
    # Functions requiring 1 parameter
    if ("mWeight" in given_parameters):
        results["BMR_Kleiber"] = BMR_Kleiber(mWeight=measurements["mWeight"])
        BMR.append(results["BMR_Kleiber"])
    if ("mHeight" in given_parameters):
        results["IBW_Lemmens"] = IBW_Lemmens(mHeight=measurements["mHeight"])
        IBW.append(results["IBW_Lemmens"])
    # Functions requiring 2 parameters
    if ("mHeight" and "cNeck") in given_parameters:
        results["AT_NHtR"] = AT_NHtR(cNeck=measurements["cNeck"], mHeight=measurements["mHeight"])
    if ("mHeight" and "cCalf") in given_parameters:
        results["AT_CHtR"] = AT_CHtR(cCalf=measurements["cCalf"], mHeight=measurements["mHeight"])
    if ("cThigh" and "cCalf") in given_parameters:
        results["AT_TCR"] = AT_TCR(cThigh=measurements["cThigh"], cCalf=measurements["cCalf"])
    if ("cForearm" and "cWrist") in given_parameters:
        results["AT_FWR"] = AT_FWR(cForearm=measurements["cForearm"], cWrist=measurements["cWrist"])
    if ("cWaist" and "cHip") in given_parameters:
        results["AI_WHR"] = AI_WHR(cWaist=measurements["cWaist"], cHip=measurements["cHip"])
    if ("cWaist" and "cCalf") in given_parameters:
        results["AI_WCR"] = AI_WCR(cWaist=measurements["cWaist"], cCalf=measurements["cCalf"])
    if ("cWaist" and "cThigh") in given_parameters:
        results["AI_WTR"] = AI_WTR(cWaist=measurements["cWaist"], cThigh=measurements["cThigh"])
    if ("mWeight" and "mHeight") in given_parameters:
        results["AT_BMI"] = AT_BMI(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"])
        results["AT_PI"] = AT_PI(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"])
        results["BSA_DuBois"] = BSA_DuBois(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"])
        results["BSA_Mosteller"] = BSA_Mosteller(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"])
        results["BSA_Haycock"] = BSA_Haycock(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"])
        results["BSA_GehanGeorge"] = BSA_GehanGeorge(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"])
        results["BSA_Fujimoto"] = BSA_Fujimoto(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"])
        BSA.append(results["BSA_DuBois"])
        BSA.append(results["BSA_Mosteller"])
        BSA.append(results["BSA_Haycock"])
        BSA.append(results["BSA_GehanGeorge"])
        BSA.append(results["BSA_Fujimoto"])
        results["IBW_Peterson"] = IBW_Peterson(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"])
        IBW.append(results["IBW_Peterson"])
    if ("mHeight" and "aGender") in given_parameters:
        results["IBW_Robinson"] = IBW_Robinson(mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        results["IBW_Miller"] = IBW_Miller(mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        results["IBW_Hamwi"] = IBW_Hamwi(mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        results["IBW_Devine"] = IBW_Devine(mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        results["IBW_KeysBrozek"] = IBW_KeysBrozek(mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        IBW.append(results["IBW_Robinson"])
        IBW.append(results["IBW_Miller"])
        IBW.append(results["IBW_Hamwi"])
        IBW.append(results["IBW_Devine"])
        IBW.append(results["IBW_KeysBrozek"])
    if ("mHeight" and "cWaist") in given_parameters:
        results["AI_WHtR"] = AI_WHtR(mHeight=measurements["mHeight"], cWaist=measurements["cWaist"])
        results["AI_WHtPI"] = AI_WHtPI(mHeight=measurements["mHeight"], cWaist=measurements["cWaist"])
        results["AI_BRI"] = AI_BRI(mHeight=measurements["mHeight"], cWaist=measurements["cWaist"])
    if ("mWeight" and "cWaist") in given_parameters:
        results["BFc_BehnkeWilmore"] = BFc_BehnkeWilmore(mWeight=measurements["mWeight"], cWaist=measurements["cWaist"])
        BF.append(results["BFc_BehnkeWilmore"])
    if ("mHeight" and "cHip") in given_parameters:
        results["BFc_BAI"] = BFc_BAI(mHeight=measurements["mHeight"], cHip=measurements["cHip"])
        BF.append(results["BFc_BAI"])
    if ("mWeight" and "sAbdominal") in given_parameters:
        results["BFs_BehnkeWilmore"] = BFs_BehnkeWilmore(mWeight=measurements["mWeight"], sAbdominal=measurements["sAbdominal"])
        BF.append(results["BFs_BehnkeWilmore"])
    # Functions requiring 3 parameters
    if ("mWeight" and "mHeight" and "cWaist") in given_parameters:
        results["AI_ABSI"] = AI_ABSI(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], cWaist=measurements["cWaist"])
        results["AI_ConicityIndex"] = AI_ConicityIndex(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], cWaist=measurements["cWaist"])
    if ("mHeight" and "cWrist" and "aGender") in given_parameters:
        results["FC_Metropolitan"] = FC_Metropolitan(mHeight=measurements["mHeight"], cWrist=measurements["cWrist"], aGender=measurements["aGender"])
    if ("mWeight" and "mHeight" and "aGender") in given_parameters:
        results["LBM_Boer"] = LBM_Boer(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        results["LBM_James"] = LBM_James(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        results["LBM_Hume"] = LBM_Hume(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        LBM.append(results["LBM_Boer"])
        LBM.append(results["LBM_James"])
        LBM.append(results["LBM_Hume"])
        results["BMR_Cunningham"] = BMR_Cunningham(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        results["BMR_KatchMcArdle"] = BMR_KatchMcArdle(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], aGender=measurements["aGender"])
        BMR.append(results["BMR_Cunningham"])
        BMR.append(results["BMR_KatchMcArdle"])
    if ("mWeight" and "cWaist" and "aGender") in given_parameters:
        results["BFc_YMCA"] = BFc_YMCA(mWeight=measurements["mWeight"], cWaist=measurements["cWaist"], aGender=measurements["aGender"])
        BF.append(results["BFc_YMCA"])
    if ("mHeight" and "cWaist" and "aGender") in given_parameters:
        results["BFc_RFM"] = BFc_RFM(mHeight=measurements["mHeight"], cWaist=measurements["cWaist"], aGender=measurements["aGender"])
        BF.append(results["BFc_RFM"])
    if ("sTricep" and "sSuprailiac" and "aGender") in given_parameters:
        results["BFs_Sloan"] = BFs_Sloan(sTricep=measurements["sTricep"], sSuprailiac=measurements["sSuprailiac"], aGender=measurements["aGender"])
        BF.append(results["BFs_Sloan"])
    if ("mWeight" and "aAge" and "aGender") in given_parameters:
        results["BMR_Schofield"] = BMR_Schofield(mWeight=measurements["mWeight"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        results["BMR_FAO2004"] = BMR_FAO2004(mWeight=measurements["mWeight"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        BMR.append(results["BMR_Schofield"])
        BMR.append(results["BMR_FAO2004"])
    # Functions requiring 4 parameters
    if ("mWeight" and "mHeight" and "aAge" and "aGender") in given_parameters:
        results["BMR_MifflinStJeor"] = BMR_MifflinStJeor(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        results["BMR_HarrisBenedict"] = BMR_HarrisBenedict(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        results["BMR_RozaShizgal"] = BMR_RozaShizgal(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        BMR.append(results["BMR_MifflinStJeor"])
        BMR.append(results["BMR_HarrisBenedict"])
        BMR.append(results["BMR_RozaShizgal"])
        results["BFx_BMI"] = BFx_BMI(aBMI=AT_BMI(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"]), aAge=measurements["aAge"], aGender=measurements["aGender"])
        BF.append(results["BFx_BMI"])
    # Functions requiring 5 parameters
    if ("mHeight" and "cNeck" and "cWaist" and "cHip" and "aGender") in given_parameters:
        results["BFc_USNavy"] = BFc_USNavy(mHeight=measurements["mHeight"], cNeck=measurements["cNeck"], cWaist=measurements["cWaist"], cHip=measurements["cHip"], aGender=measurements["aGender"])
        BF.append(results["BFc_USNavy"])
    if ("sChest" and "sAbdominal" and "sThigh" and "aAge" and "aGender") in given_parameters:
        results["BFs_JacksonPollock3"] = BFs_JacksonPollock3(sChest=measurements["sChest"], sAbdominal=measurements["sAbdominal"], sThigh=measurements["sThigh"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        BF.append(results["BFs_JacksonPollock3"])
    # Functions requiring 6 parameters
    if ("cWaist" and "mWeight" and "cWaist" and "cHip" and "cForearm" and "aGender") in measurements:
        results["BFs_JacksonPollock4"] = BFs_JacksonPollock4(sAbdominal=measurements["sAbdominal"], sTricep=measurements["sTricep"], sThigh=measurements["sThigh"], sSuprailiac=measurements["sSuprailiac"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        BF.append(results["BFs_JacksonPollock4"])
    if ("sBicep" and "sTricep" and "sSubscapular" and "sSuprailiac" and "aAge" and "aGender") in measurements:
        results["BFs_DurninWomersley"] = BFs_DurninWomersley(sBicep=measurements["sBicep"], sTricep=measurements["sTricep"], sSubscapular=measurements["sSubscapular"], sSuprailiac=measurements["sSuprailiac"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        BF.append(results["BFs_DurninWomersley"])
    if ("aAge" and "mHeight" and "sBicep" and "sTricep" and "cWaist" and "aGender") in measurements:
        results["BFx_Henry2018"] = BFx_Henry2018(aAge=measurements["aAge"], mHeight=measurements["mHeight"], sBicep=measurements["sBicep"], sTricep=measurements["sTricep"], cWaist=measurements["cWaist"], aGender=measurements["aGender"])
        BF.append(results["BFx_Henry2018"])
    # Functions requiring 8 parameters
    if ("aAge" and "cWaist" and "cHip" and "cForearm" and "cWrist"and "cThigh" and "cCalf" and "aGender") in measurements:
        results["BFc_CovertBailey"] = BFc_CovertBailey(aAge=measurements["aAge"], cWaist=measurements["cWaist"], cHip=measurements["cHip"], cForearm=measurements["cForearm"], cWrist=measurements["cWrist"], cThigh=measurements["cThigh"], cCalf=measurements["cCalf"], aGender=measurements["aGender"])
        BF.append(results["BFc_CovertBailey"])
    # Functions requiring 9 parameters
    if ("mHeight" and "mWeight" and "sTricep" and "sSubscapular" and "sSuprailiac" and "bHumerus" and "bFemur" and "cBicep" and "cCalf") in measurements:
        scores = ST_HeathCarter(mWeight=measurements["mWeight"], mHeight=measurements["mHeight"], sTricep=measurements["sTricep"], sSubscapular=measurements["sSubscapular"], sSuprailiac=measurements["sSuprailiac"], bHumerus=measurements["bHumerus"], bFemur=measurements["bFemur"], cBicep=measurements["cBicep"], cCalf=measurements["cCalf"])
        results["ST_HeathCarter (Endomorphy)"] = scores["Endomorphy"]
        results["ST_HeathCarter (Mesomorphy)"] = scores["Mesomorphy"]
        results["ST_HeathCarter (Ectomorphy)"] = scores["Ectomorphy"]
    # Functions requiring 10 parameters
    if ("mWeight" and "sChest" and "sAbdominal" and "sThigh" and "sTricep" and "sSubscapular" and "sSuprailiac" and "sMidaxillary" and "aAge" and "aGender") in measurements:
        results["BFs_JacksonPollock7"] = BFs_JacksonPollock7(mWeight=measurements["mWeight"], sChest=measurements["sChest"], sAbdominal=measurements["sAbdominal"], sThigh=measurements["sThigh"], sTricep=measurements["sTricep"], sSubscapular=measurements["sSubscapular"], sSuprailiac=measurements["sSuprailiac"], sMidaxillary=measurements["sMidaxillary"], aAge=measurements["aAge"], aGender=measurements["aGender"])
        BF.append(results["BFs_JacksonPollock7"])
    if ("mWeight" and "sChest" and "sAbdominal" and "sThigh" and "sTricep" and "sSubscapular" and "sSuprailiac" and "sLowerback" and "sCalf" and "sBicep") in measurements:
        results["BFs_Parrillo"] = BFs_Parrillo(mWeight=measurements["mWeight"], sChest=measurements["sChest"], sAbdominal=measurements["sAbdominal"], sThigh=measurements["sThigh"], sTricep=measurements["sTricep"], sSubscapular=measurements["sSubscapular"], sSuprailiac=measurements["sSuprailiac"], sLowerback=measurements["sLowerback"], sCalf=measurements["sCalf"], sBicep=measurements["sBicep"])
        BF.append(results["BFs_Parrillo"])
    if ("aAge" and "sChest" and "sAbdominal" and "sThigh" and "cWaist" and "cForearm" and "sTricep" and "sSuprailiac" and "cHip" and "aGender") in measurements:
        results["BFx_JacksonPollock3Girth"] = BFx_JacksonPollock3Girth(aAge=measurements["aAge"], sChest=measurements["sChest"], sAbdominal=measurements["sAbdominal"], sThigh=measurements["sThigh"], cWaist=measurements["cWaist"], cForearm=measurements["cForearm"], sTricep=measurements["sTricep"], sSuprailiac=measurements["sSuprailiac"], cHip=measurements["cHip"], aGender=measurements["aGender"])
        BF.append(results["BFx_JacksonPollock3Girth"])
    # Calculating averages
    results["BF-Average"] = _round(sum(BF) / len (BF), 3)
    results["BMR-Average"] = _round(sum(BMR) / len (BMR), 2)
    results["BSA-Average"] = _round(sum(BSA) / len (BSA), 4)
    results["LBM-Average"] = _round(sum(LBM) / len (LBM), 1)
    results["IBW-Average"] = _round(sum(IBW) / len (IBW), 1)
    # Return results
    return {k: v for k, v in results.items() if v is not None}


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1].lower() == "test":
        tester()
    else:
        print("Anthropy module loaded.")
        print("Run 'python anthropy.py test' to execute built-in tests.")