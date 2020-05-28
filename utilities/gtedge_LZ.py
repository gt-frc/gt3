from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from ImpRadiation.ImpurityRadiation import ImpRad


def calc_Lz_gtedge(T_ev_val):
    if T_ev_val < 3.0:
        T_ev_val = 3.0
    if T_ev_val > 10000:
        T_ev_val = 10000

    if 3.0 <= T_ev_val < 10.0:
        # dte = 0.1
        Lz_val = 1 / (7.771844618017265E+17 +
                      59122888943831.57 * np.exp(T_ev_val) +
                      2.426539413730379E+20 * np.exp(-T_ev_val))
    if 10.0 <=T_ev_val< 100.0:
        # dte = 0.5
        Lz_val = (T_ev_val*(-2.654182242083568E-20 +
                        T_ev_val*(2.458509606418251E-22 +
                                T_ev_val*(-9.488347829679447E-25 +
                                        T_ev_val*(-1.316645200326443E-27 +
                                                T_ev_val*(1.459022440553818E-29))))) +
                 1.639722967154732E-18 +
                  (1/T_ev_val)*(-6.374237210938420E-17 +
                        (1/T_ev_val)*(1.590289149116386E-15 +
                               (1/T_ev_val)*(-2.435345151607877E-14 +
                                      (1/T_ev_val)*(2.038708041891294E-13 +
                                              (1 / T_ev_val)*(-6.468093918696965E-13))))))
    if 100.0 <=T_ev_val< 10000.0:
        # if 100.0 <=T_ev_val< 120.0:
        #     dte = 0.5
        # if 120.0 <=T_ev_val< 1000.0:
        #     dte = 5.0
        # if 1000.0 <=T_ev_val< 10000.0:
        #     dte = 10.0

        Lz_val = (T_ev_val*(5.333211787469455E-25+
                      T_ev_val*(-1.404553291126376E-28+
                         T_ev_val*(2.164576220631324E-32+
                            T_ev_val*(-1.691926173232754E-36+
                               T_ev_val*(5.202445752080436E-41)))))
                - 5.770913821900679E-22 +
                (1/T_ev_val) * (1.042550985018631E-18+
                     (1/T_ev_val)*(-3.134926132115641E-16+
                        (1/T_ev_val)*(8.184722721684790E-14+
                           (1/T_ev_val)*(-6.157523360890277E-12+
                              (1/T_ev_val)*(9.848433223044370E-11))))))

    Te_ev_min = np.log10(3.0)
    Te_ev_max = np.log10(10000.0)
    Te_ev = np.logspace(Te_ev_min, Te_ev_max, 1000)

    Lz = np.zeros(np.asarray(Te_ev).shape)
    for i,T in enumerate(np.asarray(Te_ev)):
        if T < 3.0:
            T = 3.0
        if T > 10000:
            T = 10000

        if 3.0 <= T < 10.0:
            # dte = 0.1
            Lz[i] = 1 / (7.771844618017265E+17 +
                          59122888943831.57 * np.exp(T) +
                          2.426539413730379E+20 * np.exp(-T))
        if 10.0 <= T < 100.0:
            # dte = 0.5
            Lz[i] = (T*(-2.654182242083568E-20 +
                        T*(2.458509606418251E-22 +
                           T*(-9.488347829679447E-25 +
                              T*(-1.316645200326443E-27 +
                                 T*(1.459022440553818E-29))))) +
                      1.639722967154732E-18 +
                      (1/T)*(-6.374237210938420E-17 +
                             (1/T)*(1.590289149116386E-15 +
                                    (1/T)*(-2.435345151607877E-14 +
                                           (1/T)*(2.038708041891294E-13 +
                                                  (1/T)*(-6.468093918696965E-13))))))
        if 100.0 <= T < 10000.0:
            # if 100.0 <=T< 120.0:
            #     dte = 0.5
            # if 120.0 <=T< 1000.0:
            #     dte = 5.0
            # if 1000.0 <=T< 10000.0:
            #     dte = 10.0

            Lz[i] = (T * (5.333211787469455E-25 +
                          T * (-1.404553291126376E-28 +
                               T * (2.164576220631324E-32 +
                                    T * (-1.691926173232754E-36 +
                                         T * (5.202445752080436E-41)))))
                      - 5.770913821900679E-22 +
                      (1/T)*(1.042550985018631E-18 +
                             (1/T)*(-3.134926132115641E-16 +
                                    (1/T)*(8.184722721684790E-14 +
                                           (1/T)*(-6.157523360890277E-12 +
                                                  (1/T)*(9.848433223044370E-11))))))

    #plt.loglog(Te_ev, UnivariateSpline(Te_ev*1.6021E-19, Lz, k=1, s=0).derivative()(Te_ev*1.6021E-19))
    #plt.show()
    dLzdT_val = UnivariateSpline(Te_ev*1.6021E-19, Lz, k=1, s=0).derivative()(T_ev_val*1.6021E-19)
    print
    print 'T_ev_val = ',T_ev_val
    print 'Lz_val = ',Lz_val
    print 'dLzdT_val = ',dLzdT_val

    print
    return Lz_val, dLzdT_val

def calc_Lz_gtedge_arr(Te_ev):

    Lz = np.zeros(np.asarray(Te_ev).shape)
    for i,T in enumerate(np.asarray(Te_ev)):
        if T < 3.0:
            T = 3.0
        if T > 10000:
            T = 10000

        if 3.0 <= T < 10.0:
            # dte = 0.1
            Lz[i] = 1 / (7.771844618017265E+17 + 59122888943831.57 * np.exp(T) + 2.426539413730379E+20 * np.exp(-T))
        if 10.0 <= T < 100.0:
            # dte = 0.5
            Lz[i] = T * (-2.654182242083568E-20 + T * (2.458509606418251E-22 + T * (
                -9.488347829679447E-25+T * (-1.316645200326443E-27+T * (
                1.459022440553818E-29))))) + 1.639722967154732E-18 + (1.0 / T) * (
                -6.374237210938420E-17+(1.0 / T) * (1.590289149116386E-15+(1.0 / T) * (-
                2.435345151607877E-14+(1.0 / T) * (2.038708041891294E-13+(1.0 / T) * (-
                6.468093918696965E-13)))))
        if 100.0 <= T < 10000.0:
            # if 100.0 <= T < 120.0:
            #     dte = 0.5
            # if 120.0 <= T < 1000.0:
            #     dte = 5.0
            # if 1000.0 <= T < 10000.0:
            #     dte = 10.0

            Lz[i] = T * (5.333211787469455E-25 + T * (-1.404553291126376E-28 + T * (
                2.164576220631324E-32+T * (-1.691926173232754E-36+T * (
                5.202445752080436E-41)))))  -5.770913821900679E-22 + (1.0 / T) * (
                1.042550985018631E-18+(1.0 / T) * (-3.134926132115641E-16+(1.0 / T) * (
                8.184722721684790E-14+(1.0 / T) * (-6.157523360890277E-12+(1.0 / T) * (
                9.848433223044370E-11)))))
    dLzdT = UnivariateSpline(Te_ev*1.6021E-19, Lz, k=1, s=0).derivative()(Te_ev*1.6021E-19)
    return Lz, dLzdT


if __name__ == '__main__':
    Te_ev_min = np.log10(3.0)
    Te_ev_max = np.log10(10000.0)
    Te_ev = np.logspace(Te_ev_min, Te_ev_max,1000)
    Lz, dLzdT = calc_Lz_gtedge_arr(Te_ev)
    plt.loglog(Te_ev, Lz*1E-13, label='gtedge')
    #plt.loglog(Te_ev, dLzdT*1E-13, label='gtedge')

    C_6 = ImpRad(z=6)
    Te_kev = Te_ev *1E-3
    Lz_gt3 = C_6.Lz(np.log10(0.002), np.log10(1E-7), np.log10(Te_kev))
    dLzdT_gt3 = C_6.dLzdT(np.log10(0.002), np.log10(1E-7), np.log10(Te_kev))
    plt.loglog(Te_ev, Lz_gt3,label='gt3')
    #plt.loglog(Te_ev, dLzdT_gt3,label='gt3')
    plt.legend()
    plt.show()