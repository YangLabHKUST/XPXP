import pandas as pd
import numpy as np
import scipy.stats as st
from scipy.stats import norm
import argparse

def calculate_p(z_scores: np.array) -> np.array:
    # Calculate constants used in determination of P values
    ln = np.log  # pylint: disable=invalid-name
    LN_2 = ln(2.0)
    RECIP_LN_10 = np.reciprocal(ln(10.0))
    """
    Function that calculates P for the TRAM results
    :param z_scores: Z scores
    :return: P values 
             (as strings, to allow for very large negative exponents)
    """
    # Since P = 2 * normal_cdf(-|Z|), P = e ^ (log_normal_cdf(-|Z|) + ln 2)
    # This can be changed to base 10 as P = 10 ^ ((log_normal_cdf(-|Z|) + ln 2) / ln 10)
    log_10_p = RECIP_LN_10 * (norm.logcdf(-np.abs(z_scores)) + LN_2)

    # Break up the log based 10 of P values into the integer and fractional part
    # To handle the case of Z = 0 (and not result in "10e-1"), set initial values to (-1.0, 1.0)
    frac_part, int_part = np.full_like(z_scores, -1.0), np.full_like(z_scores, 1.0)
    np.modf(log_10_p, out=(frac_part, int_part), where=(z_scores != 0.0))

    # Construct strings for the P values
    # 1) Add one to the fractional part to ensure that the result mantissa is between 1 and 10
    # 2) Subtract one from the integer part to compensate and keep the overall value correct
    result = np.char.add(np.char.add(np.power(10.0, (frac_part + 1.0)).astype(str), 'e'),
                         (int_part - 1).astype(int).astype(str))

    return result


if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='reformat ss')
    parser.add_argument('--inp', type=str, help='output file path', required=True)
    args = parser.parse_args()
    
    df = pd.read_csv(args.inp,sep='\t').dropna()
    df['P'] = calculate_p(df['Z'].values).astype(float)
    df['BETA'] = df['Z']
    df[['SNP','A1','A2','BETA','P']].to_csv(args.inp+'.prscsx_format.txt',sep='\t',index=None)  