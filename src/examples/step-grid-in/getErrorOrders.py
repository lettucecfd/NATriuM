import glob

directory = "/work/ws-tmp/ps500501-natrium/"
# directory = "/mnt/c/Users/phili/Desktop/"
directories = glob.glob(directory + "*/*.out")
directories.sort()

for filename in directories:
    runname = filename.split("/")[-2]
    with open(filename) as file:
        string = file.read().splitlines()
        for line in string:
            if "Actual dt:" in line:
                dt = float(line.removeprefix("::::Actual dt:                ").removesuffix(" s"))
            elif "dx_min:" in line:
                dxmin = float(line.removeprefix("::::dx_min:                   "))
            elif "dx_max:" in line:
                dxmax = float(line.removeprefix("::::dx_max:                   "))
        p = 4
        dx = dxmin
        Omin = min(dx**(p+1)/dt,dx**p)+dt**2
        dx = dxmax
        Omax = min(dx**(p+1)/dt,dx**p)+dt**2
        output = f"{runname}: Omin={Omin:.3e},\tOmax={Omax:.3e}\t"
        if dxmin**p < dxmin**(p+1)/dxmin:
            output += "dxmin dominates"
        else:
            output += "dt dominates dxmin"
        output += "\t"
        if dxmax**p < dxmax**(p+1)/dxmax:
            output += "dxmax dominates"
        else:
            output += "dt dominates dxmax"
        print(output)



string = """9671280_cfl4_ref0_glc_p4_newNonUni: Omin=1.189e-08,     Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9677567_cfl4_ref0_glc_p4_newNonUni: Omin=8.524e-09,     Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9682930_cfl4_ref__p_newNonUni: Omin=8.494e-09,  Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9684155_cfl4_ref__p_newNonUni: Omin=8.494e-09,  Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9684384_cfl4_ref__p_newUni: Omin=8.494e-09,     Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9684432_cfl2_ref__p_newNonUni: Omin=2.124e-09,  Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9684433_cfl3_ref__p_newNonUni: Omin=4.778e-09,  Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9684434_cfl5_ref__p_newNonUni: Omin=1.327e-08,  Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9684476_cfl4_ref__p_cartesianNonUni: Omin=1.327e-08,    Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9684651_cfl4_ref__p_newUni: Omin=3.874e-08,     Omax=1.507e-02  dt dominates dxmin      dt dominates dxmax
9684661_cfl1.2_ref__p_newUniCoarsen3: Omin=2.361e-06,   Omax=2.562e-01  dt dominates dxmin      dxmax dominates
9684663_cfl4_ref__p_newNonUniCoarse2: Omin=3.398e-08,   Omax=1.494e-01  dt dominates dxmin      dt dominates dxmax
9684664_cfl4_ref__p_newNonUni: Omin=2.123e-09,  Omax=5.836e-04  dt dominates dxmin      dt dominates dxmax
9684665_cfl2_ref__p_newNonUniCoarse2: Omin=8.495e-09,   Omax=1.494e-01  dt dominates dxmin      dt dominates dxmax
9686957_cfl4_ref1__p_newNonUni: Omin=2.123e-09, Omax=5.836e-04  dt dominates dxmin      dt dominates dxmax
9687124_cfl4_ref0__p5_newNonUni_0: Omin=3.479e-09,      Omax=7.572e-03  dt dominates dxmin      dt dominates dxmax
9687198_cfl4_ref0__p4_newNonUni_1: Omin=3.398e-08,      Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687199_cfl4_ref0__p4_newNonUni_2: Omin=3.398e-08,      Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687200_cfl4_ref0__p4_newNonUni_3: Omin=3.398e-08,      Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687201_cfl4_ref0__p4_newNonUni_4: Omin=3.398e-08,      Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687206_cfl4_ref0__p6_newNonUni_0: Omin=6.712e-09,      Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687720_cfl4_ref0__p4_nonUniRes10_0: Omin=6.712e-09,    Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687721_cfl4_ref0__p4_nonUniRes20_0: Omin=6.712e-09,    Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687722_cfl4_ref0__p4_nonUniRes25_0: Omin=6.712e-09,    Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687723_cfl4_ref0__p4_nonUniRes30_0: Omin=6.712e-09,    Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9687724_cfl4_ref0__p4_nonUniRes40_0: Omin=7.957e-08,    Omax=2.519e-02  dt dominates dxmin      dt dominates dxmax
9687725_cfl4_ref0__p4_nonUniRes50_0: Omin=3.264e-08,    Omax=1.359e-01  dt dominates dxmin      dxmax dominates
9687726_cfl4_ref0__p4_nonUniRes60_0: Omin=2.514e-07,    Omax=2.574e-02  dt dominates dxmin      dt dominates dxmax
9687727_cfl4_ref0__p4_nonUniRes75_0: Omin=1.031e-07,    Omax=2.073e-02  dt dominates dxmin      dt dominates dxmax
9687728_cfl4_ref0__p4_nonUniRes80_0: Omin=7.957e-08,    Omax=2.101e-02  dt dominates dxmin      dt dominates dxmax
9687729_cfl4_ref0__p4_nonUniRes100_0: Omin=3.264e-08,   Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9689791_cfl4_ref0__p4_newNonUni_0deg_Re20000: Omin=3.398e-08,   Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9689792_cfl4_ref0__p4_newNonUni_0deg_Re50000: Omin=3.398e-08,   Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9689793_cfl4_ref0__p4_newNonUni_0deg_Re100000: Omin=3.398e-08,  Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9690806_cfl4_ref0__p4_nonUniRes100_10deg: Omin=3.262e-08,       Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9693817_cfl4_ref0__p4_nonUniRes100_0deg_Ma0.7: Omin=7.108e-09,  Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9693818_cfl4_ref0__p4_nonUniRes100_0deg_Ma0.9: Omin=1.175e-08,  Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9693848_cfl4_ref1_newNonUni_0deg_Re20000: Omin=8.494e-09,       Omax=1.908e-03  dt dominates dxmin      dt dominates dxmax
9693849_cfl4_ref1_newNonUni_0deg_Re50000: Omin=8.494e-09,       Omax=1.908e-03  dt dominates dxmin      dt dominates dxmax
9694158_cfl4_ref0_newNonUniRes100_0deg_Re10000: Omin=2.555e-07, Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9694159_cfl4_ref0_newNonUniRes100_10deg_Re10000: Omin=2.516e-07,        Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9694162_cfl4_ref0_newNonUniRes100_20deg_Re10000: Omin=2.403e-07,        Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9695214_cfl4_ref0_03_newCartesian_0deg_Re10000: Omin=7.463e-08, Omax=1.393e-04  dt dominates dxmin      dt dominates dxmax
9695341_cfl3_ref1_03_newCartesian_0deg_Re10000: Omin=1.049e-08, Omax=8.714e-06  dt dominates dxmin      dt dominates dxmax
9695342_cfl4_ref1_03_newCartesian_0deg_Re10000: Omin=1.866e-08, Omax=8.722e-06  dt dominates dxmin      dt dominates dxmax
9695346_cfl4_ref0_nonUniRectangular_0deg_Re10000: Omin=2.555e-07,       Omax=1.096e-05  dt dominates dxmin      dt dominates dxmax
9695347_cfl4_ref1_nonUniRectangular_0deg_Re10000: Omin=6.386e-08,       Omax=7.337e-07  dt dominates dxmin      dxmax dominates
9695357_cfl1_ref0_03_newCartesian_0deg_Re10000: Omin=4.669e-09, Omax=1.393e-04  dt dominates dxmin      dt dominates dxmax
9695358_cfl2_ref0_03_newCartesian_0deg_Re10000: Omin=1.866e-08, Omax=1.393e-04  dt dominates dxmin      dt dominates dxmax
9695360_cfl6_ref0_03_newCartesian_0deg_Re10000: Omin=1.679e-07, Omax=1.394e-04  dt dominates dxmin      dt dominates dxmax
9695361_cfl8_ref0_03_newCartesian_0deg_Re10000: Omin=2.985e-07, Omax=1.396e-04  dt dominates dxmin      dt dominates dxmax
9695387_cfl1_ref0_p4_nonUniRectangular: Omin=1.602e-08, Omax=1.072e-05  dt dominates dxmin      dt dominates dxmax
9695388_cfl1_ref0_p5_nonUniRectangular: Omin=6.599e-09, Omax=1.071e-05  dt dominates dxmin      dt dominates dxmax
9695389_cfl1_ref1_p4_nonUniRectangular: Omin=3.995e-09, Omax=6.739e-07  dt dominates dxmin      dxmax dominates
9695390_cfl1_ref1_p5_nonUniRectangular: Omin=1.638e-09, Omax=6.715e-07  dt dominates dxmin      dxmax dominates
9695392_cfl1_ref2_p5_nonUniRectangular: Omin=1.638e-09, Omax=6.715e-07  dt dominates dxmin      dxmax dominates
9695393_cfl4_ref0_p4_nonUniRectangular: Omin=2.555e-07, Omax=1.096e-05  dt dominates dxmin      dt dominates dxmax
9695394_cfl4_ref0_p5_nonUniRectangular: Omin=1.047e-07, Omax=1.081e-05  dt dominates dxmin      dt dominates dxmax
9695396_cfl4_ref1_p5_nonUniRectangular: Omin=2.616e-08, Omax=6.960e-07  dt dominates dxmin      dxmax dominates
9695398_cfl4_ref2_p5_nonUniRectangular: Omin=2.616e-08, Omax=6.960e-07  dt dominates dxmin      dxmax dominates
9695399_cfl8_ref0_p4_nonUniRectangular: Omin=1.022e-06, Omax=1.173e-05  dt dominates dxmin      dt dominates dxmax
9695400_cfl8_ref0_p5_nonUniRectangular: Omin=4.185e-07, Omax=1.112e-05  dt dominates dxmin      dt dominates dxmax
9695401_cfl8_ref1_p4_nonUniRectangular: Omin=2.554e-07, Omax=9.253e-07  dt dominates dxmin      dxmax dominates
9695402_cfl8_ref1_p5_nonUniRectangular: Omin=1.046e-07, Omax=7.745e-07  dt dominates dxmin      dxmax dominates
9695403_cfl8_ref2_p4_nonUniRectangular: Omin=6.385e-08, Omax=1.058e-07  dt dominates dxmin      dxmax dominates
9695404_cfl8_ref2_p5_nonUniRectangular: Omin=6.385e-08, Omax=1.058e-07  dt dominates dxmin      dxmax dominates
9695525_cfl4_ref0_03_newCartesianRegular_0deg_Re10000: Omin=3.291e-06,  Omax=3.295e-06  dxmin dominates dt dominates dxmax
9695526_cfl4_ref1_03_newCartesianRegular_0deg_Re10000: Omin=8.209e-07,  Omax=8.212e-07  dxmin dominates dt dominates dxmax
9695534_cfl4_ref3_03_newCartesianRegular_0deg_Re10000: Omin=8.209e-07,  Omax=8.212e-07  dxmin dominates dt dominates dxmax
9695535_cfl4_ref4_03_newCartesianRegular_0deg_Re10000: Omin=8.209e-07,  Omax=8.212e-07  dxmin dominates dt dominates dxmax
9695558_cfl4_ref2_03_newCartesianRegular_0deg_Re10000: Omin=8.209e-07,  Omax=8.212e-07  dxmin dominates dt dominates dxmax
9695560_cfl4_ref0_nonUniRes100_Ma0.5: Omin=3.627e-09,   Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9695561_cfl4_ref0_nonUniRes100_Ma0.7: Omin=7.108e-09,   Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9695562_cfl4_ref0_nonUniRes100_Ma0.9: Omin=1.175e-08,   Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9695571_cfl2_ref0_p4_nonUniRectangular: Omin=6.391e-08, Omax=1.077e-05  dt dominates dxmin      dt dominates dxmax
9695574_cfl5_ref0_p4_nonUniRectangular: Omin=3.991e-07, Omax=4.622e-04  dt dominates dxmin      dt dominates dxmax
9696188_cfl3_ref0_p4_nonUniRectangular: Omin=1.437e-07, Omax=4.619e-04  dt dominates dxmin      dt dominates dxmax
9696191_cfl4_ref0_p4_nonUniRectangularStraightInlet: Omin=3.185e-08,    Omax=5.117e-04  dt dominates dxmin      dt dominates dxmax
9696451_cfl4_ref0_newNonUni_0deg_Re10000: Omin=3.058e-08,       Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9696452_cfl4_ref0_nonUniRectangular_0deg_Re10000: Omin=2.299e-07,       Omax=4.620e-04  dt dominates dxmin      dt dominates dxmax
9696453_cfl4_ref0_nonUniRectangularStraightInlet_0deg_Re10000: Omin=2.866e-08,  Omax=5.117e-04  dt dominates dxmin      dt dominates dxmax
9696468_cfl4_ref0_nonUniRectangularStraightInletZGWall_0deg_Re10000: Omin=2.299e-07,    Omax=5.928e-04  dt dominates dxmin      dt dominates dxmax
9696523_cfl4_ref0_nonUniRectangular_0deg_Re10000: Omin=2.299e-07,       Omax=4.620e-04  dt dominates dxmin      dt dominates dxmax
9696619_cfl4_ref0_nonUniRectangularStraightInlet_0deg_Re10000: Omin=2.866e-08,  Omax=5.117e-04  dt dominates dxmin      dt dominates dxmax
9696661_cfl4_ref0_nonUniRectangularStraightInletZGWall_0deg_Re10000: Omin=2.299e-07,    Omax=5.928e-04  dt dominates dxmin      dt dominates dxmax
9697452_cfl4_ref0_nonUniRectangularZGWall_0deg_Re10000: Omin=2.299e-07, Omax=4.620e-04  dt dominates dxmin      dt dominates dxmax
9697454_cfl4_ref1_newNonUniTrip_0deg_Re10000: Omin=1.462e-08,   Omax=4.800e-03  dt dominates dxmin      dt dominates dxmax
9697473_cfl4_ref0_newNonUniTrip_0deg_Re10000: Omin=3.736e-08,   Omax=4.683e-02  dt dominates dxmin      dt dominates dxmax
9697692_cfl4_ref1_newNonUniTrip2_0deg_Re10000: Omin=9.341e-09,  Omax=4.800e-03  dt dominates dxmin      dt dominates dxmax
9697693_cfl4_ref2_newNonUniTrip_0deg_Re10000: Omin=2.335e-09,   Omax=4.554e-04  dt dominates dxmin      dxmax dominates
9697708_cfl4_ref0_newNonUni_0deg_Re10000: Omin=1.359e-08,       Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9697709_cfl4_ref0_newNonUni_0deg_Re10000: Omin=1.957e-08,       Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9697713_cfl4_ref0_newUniStepped_0deg_Re10000: Omin=1.957e-08,   Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9697725_cfl4_ref0_newNonUni_0deg_Re5000: Omin=3.058e-08,        Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9697873_cfl4_ref0_nonUniRes100_0deg_Re5000: Omin=3.264e-08,     Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9697878_cfl3_ref0_nonUniRes100_0deg_Re5000: Omin=1.836e-08,     Omax=1.439e-02  dt dominates dxmin      dt dominates dxmax
9697913_cfl5_ref0_newNonUni_0deg_Re5000: Omin=5.309e-08,        Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9697914_cfl6_ref0_newNonUni_0deg_Re5000: Omin=7.645e-08,        Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9697915_cfl4_ref0_newNonUni_0deg_Re5000: Omin=6.040e-08,        Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9697916_cfl4_ref1_newNonUni_0deg_Re10000: Omin=2.416e-09,       Omax=1.908e-03  dt dominates dxmin      dt dominates dxmax
9700609_cfl10_ref0_newNonUni_0deg_Re10000: Omin=2.123e-07,      Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9700697_cfl5_ref0_newUniStepped_0deg_Re10000_Ma1.5: Omin=2.123e-07,     Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9700699_cfl5_ref1_newUniStepped_0deg_Re10000_Ma1.5: Omin=2.997e-09,     Omax=5.992e-03  dt dominates dxmin      dt dominates dxmax
9700700_cfl10_ref1_newUniStepped_0deg_Re10000_Ma1.5: Omin=2.997e-09,    Omax=5.992e-03  dt dominates dxmin      dt dominates dxmax
9700706_cfl2_ref1_newUniStepped_0deg_Re10000_Ma1.5: Omin=4.796e-10,     Omax=5.992e-03  dt dominates dxmin      dt dominates dxmax
9700716_cfl1_ref1_newUniStepped_0deg_Re10000_Ma1.5: Omin=1.199e-10,     Omax=5.992e-03  dt dominates dxmin      dt dominates dxmax
9700720_cfl20_ref0_newNonUni_0deg_Re10000_Ma1.5: Omin=1.199e-10,        Omax=5.992e-03  dt dominates dxmin      dt dominates dxmax
9700721_cfl40_ref0_newNonUni_0deg_Re10000_Ma1.5: Omin=1.199e-10,        Omax=5.992e-03  dt dominates dxmin      dt dominates dxmax
9700727_cfl10_ref1_newNonUni_0deg_Re10000_Ma1.5: Omin=5.309e-08,        Omax=1.908e-03  dt dominates dxmin      dt dominates dxmax
9702207_cfl4_ref0_newNonUni_0deg_Re300_Ma1.5: Omin=3.398e-08,   Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9702340_cfl20_ref1_newNonUni_0deg_Re300_Ma1.5: Omin=3.398e-08,  Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9702342_cfl10_ref1_sphereTransfinite2_0deg_Re10000_Ma1.5: Omin=2.319e-08,       Omax=7.605e-03  dt dominates dxmin      dt dominates dxmax
9702343_cfl20_ref1_sphereTransfinite2_0deg_Re10000_Ma1.5: Omin=2.319e-08,       Omax=7.605e-03  dt dominates dxmin      dt dominates dxmax
9702346_cfl20_ref0_newNonUni_0deg_Re300_Ma1.5: Omin=2.319e-08,  Omax=7.605e-03  dt dominates dxmin      dt dominates dxmax
9702347_cfl1_ref1_sphereTransfinite2_0deg_Re10000_Ma1.5: Omin=2.319e-10,        Omax=7.605e-03  dt dominates dxmin      dt dominates dxmax
9702348_cfl4_ref2_sphereTransfinite2_0deg_Re10000_Ma1.5: Omin=9.277e-10,        Omax=4.860e-04  dt dominates dxmin      dxmax dominates
9702349_cfl4_ref0_newNonUni_0deg_Re1000_Ma1.5: Omin=3.398e-08,  Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9702351_cfl4_ref0_newNonUni_0deg_Re5000_Ma1.5: Omin=3.398e-08,  Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9703696_cfl1_ref0_newNonUni_0deg_Re1000_Ma1.5: Omin=2.125e-09,  Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9703697_cfl0.1_ref0_newNonUni_0deg_Re1000_Ma1.5: Omin=2.231e-11,        Omax=1.670e-02  dt dominates dxmin      dt dominates dxmax
9703701_cfl4_ref1_sphereTransfinite2_0deg_Re10000_Ma1.5: Omin=3.327e-08,        Omax=7.605e-03  dt dominates dxmin      dt dominates dxmax
9704184_cfl4_ref1_sphereTransfinite2_0deg_Re10000_Ma1.5: Omin=3.327e-08,        Omax=7.605e-03  dt dominates dxmin      dt dominates dxmax
9704185_cfl4_ref2_sphereTransfinite2_0deg_Re10000_Ma1.5: Omin=8.317e-09,        Omax=4.860e-04  dt dominates dxmin      dxmax dominates"""

string = string.splitlines()
for line in string:
    job = line.split(" ")[0]
    out = job + "\n" + "Omin"
    for l in line.split("Omin")[1:]:
        out += l.replace("  ", " ").replace("  ", " ").replace("  ", " ").replace(" ","\t").replace("dt\t","dt ").replace("dominates\t","dominates ")
    print(out)