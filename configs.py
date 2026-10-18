#!/usr/bin/env python3
"""
configs.py -- constants for each of the 14 eta notebooks.

All values were extracted directly from the .nb files, not retyped.

Region naming is the notebook 2-14 convention (inner / middle / outer).
ETA 1 NOTE: eta___1.nb names the same three regions (x2bar, x2, x1) with
constants (Q2bar, Q2, Q1) and occupancies (n2bar g2bar, n2 g2, n1 g1), and
calls its couplings beta12bar / gamma12. The mapping used below is

    inner  : n2,g2,Q2 <- n2bar, g2bar, Q2bar
    middle : n3,g3,Q3 <- n2,    g2,    Q2
    outer  : n1,g1,Q1 <- n1,    g1,    Q1
    beta12 <- beta12bar        gamma12bar <- gamma12

Notebook 1 also sets aa = zz/bb = 306/1553 rather than the literal 0.197
used by notebooks 2-14; that is carried through below.
"""

from tf_core import Config

CONFIGS = {
    1: Config(name='eta___1',  n1=1,  g1=3, n2=1, g2=1, n3=1,  g3=1,
              gamma12bar=0.3954232044351274, beta12=0.2753036090258285,
              beta22bar=0.9621573728871041,
              Q1=-0.45, Q2=-2.0971474564888006, Q3=-0.851982694608053,
              m2_start=9, aa=306.0 / 1553.0),

    2: Config(name='eta____2', n1=2,  g1=3, n2=1, g2=2, n3=1,  g3=2,
              gamma12bar=0.4974508021184714, beta12=0.2485750240167111,
              beta22bar=0.882615032977455,
              Q1=-0.5, Q2=-3.6375498477951878, Q3=-1.4247837854483154,
              m2_start=14),

    3: Config(name='eta____3', n1=3,  g1=3, n2=3, g2=1, n3=3,  g3=1,
              gamma12bar=0.3954232044351274, beta12=0.22559600009756925,
              beta22bar=0.859355921116699,
              Q1=-0.5142857142857142, Q2=-2.007183585486792,
              Q3=-0.973694508123489, m2_start=22),

    4: Config(name='eta____4', n1=3,  g1=4, n2=2, g2=2, n3=2,  g3=2,
              gamma12bar=0.45566826601850596, beta12=0.2204674444696147,
              beta22bar=0.8237105803142551,
              Q1=-0.6947368421052631, Q2=-3.7583491797205992,
              Q3=-1.577809540840739, m2_start=25),

    5: Config(name='eta____5', n1=5,  g1=3, n2=5, g2=1, n3=5,  g3=1,
              gamma12bar=0.3987206213820558, beta12=0.21393789137970892,
              beta22bar=0.824988803688124,
              Q1=-0.525, Q2=-1.9948905987395598,
              Q3=-0.9585131676908107, m2_start=32),

    6: Config(name='eta____6', n1=6,  g1=3, n2=3, g2=2, n3=3,  g3=2,
              gamma12bar=0.505242178637937, beta12=0.21321181216732096,
              beta22bar=0.7946880881085071,
              Q1=-0.5275862068965518, Q2=-3.6401613735622833,
              Q3=-1.357809759306383, m2_start=30),

    7: Config(name='eta____7', n1=7,  g1=3, n2=7, g2=1, n3=7,  g3=1,
              gamma12bar=0.400144239328412, beta12=0.20902489540234132,
              beta22bar=0.81049344303531,
              Q1=-0.5294117647058824, Q2=-1.990947646406977,
              Q3=-0.9523141799457909, m2_start=41),

    8: Config(name='eta____8', n1=6,  g1=4, n2=4, g2=2, n3=4,  g3=2,
              gamma12bar=0.4586778905637988, beta12=0.20850380439721672,
              beta22bar=0.7917378610697486,
              Q1=-0.7076923076923077, Q2=-3.7573992374734764,
              Q3=-1.5526090060270206, m2_start=40),

    9: Config(name='eta____9', n1=9,  g1=3, n2=9, g2=1, n3=9,  g3=1,
              gamma12bar=0.40093784733070004, beta12=0.20631963826157407,
              beta22bar=0.8025057114461955,
              Q1=-0.5318181818181819, Q2=-1.989094374486371,
              Q3=-0.9489464140282534, m2_start=50),

    10: Config(name='eta____10', n1=10, g1=3, n2=5, g2=2, n3=5, g3=2,
               gamma12bar=0.5068195691852881, beta12=0.20660077501801427,
               beta22bar=0.7778524931404934,
               Q1=-0.5326530612244897, Q2=-3.6465672354338867,
               Q3=-1.3456132586363403, m2_start=45),

    11: Config(name='eta____11', n1=11, g1=3, n2=11, g2=1, n3=11, g3=1,
               gamma12bar=0.401443878585218, beta12=0.2046075687172968,
               beta22bar=0.7974477329126087,
               Q1=-0.5333333333333333, Q2=-1.9880392139611935,
               Q3=-0.9468309809277188, m2_start=45),

    12: Config(name='eta____12', n1=9, g1=4, n2=6, g2=2, n3=6, g3=2,
               gamma12bar=0.45968672813171524, beta12=0.20461603041128879,
               beta22bar=0.7812734132661153,
               Q1=-0.711864406779661, Q2=-3.758645495522109,
               Q3=-1.5445330058685771, m2_start=45),

    13: Config(name='eta____13', n1=13, g1=3, n2=13, g2=1, n3=13, g3=1,
               gamma12bar=0.40179466729030605, beta12=0.2034267249755678,
               beta22bar=0.7939577368574934,
               Q1=-0.534375, Q2=-1.987364953721077,
               Q3=-0.9453789011026564, m2_start=45),

    14: Config(name='eta____14', n1=14, g1=3, n2=7, g2=2, n3=7, g3=2,
               gamma12bar=0.5074975496262601, beta12=0.2038185199059122,
               beta22bar=0.770721910674057,
               Q1=-0.5347826086956522, Q2=-3.649874507155945,
               Q3=-1.34049791878099, m2_start=45),
}

if __name__ == '__main__':
    print(f'{"eta":>4} {"file":<12} {"N":>4} {"n1 g1":>7} {"n2 g2":>7} {"n3 g3":>7} {"R":>14}')
    from tf_core import Model
    for k in sorted(CONFIGS):
        c = CONFIGS[k]
        m = Model(c)
        print(f'{k:4d} {c.name:<12} {c.N:4d} {c.n1:3d}{c.g1:4d} {c.n2:3d}{c.g2:4d} '
              f'{c.n3:3d}{c.g3:4d} {m.R:14.9f}')
