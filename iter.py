import os
import shutil
import h5py
import numpy as np
import re

# Step.0, 获取当前脚本所在位置
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)


# Step.I, 读取 coupling.dat 中的迭代信息，并更新runFluent.jou
coupling_path = os.path.join(script_dir, 'coupling.dat')
with open(coupling_path, 'r') as coupling_file:
    lines = coupling_file.readlines()
    for i, line in enumerate(lines):
        if line.startswith('#Multilevel Flag'):
            flag = int(lines[i].split()[-1])
        elif line.startswith('#Omega'):
            Omega = float(lines[i].split()[-1])
        elif line.startswith('#Residual'):
            re_factor = float(lines[i].split()[-1])
        elif line.startswith('#Parallel processes for RMC'):
            processes_RMC = int(lines[i].split()[-1])
        elif line.startswith('#Parallel processes for Fluent'):
            processes_CFD = int(lines[i].split()[-1])
        elif line.startswith('#SOR flag for th_calculation'):
            SOR_flag_for_th = int(lines[i].split()[-1])
        elif line.startswith('#Output all iter flag'):
            all_iter_flag = int(lines[i].split()[-1])
        elif line.startswith('#Project Name'):
            name = lines[i+1].split()
    max_iter = int(lines[0])

jou_path = os.path.join(script_dir, 'runFluent.jou')  # runFluent.jou 文件的路径
temp_file_path = jou_path + '.tmp'  # 创建一个临时文件用于修改

with open(jou_path, 'r', encoding='utf-8') as jou_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in jou_file:
        if line.startswith('rc'):
            temp_file.write(f'rc {name[1]}.cas.h5\n')
        elif line.startswith('wd'):
            temp_file.write(f'wd {name[1]}.dat.h5\n')
        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, jou_path)  # 强制替换原始文件

os.remove(temp_file_path)  # 删除临时文件


# Step.II, 初始化命令与残差向量，并设置源文件和目标路径
run_MC = f'mpiexec -n {processes_RMC} RMC {name[0]}.rmc'
run_CFD = f'fluent 3ddp -g -i runFluent.jou -t{processes_CFD}'
rms = []
re_ave = []
keff = []
std = []
keff_re_diff = []

####################################### change this #######################################
source_files = [
    "info_fuel.h5",
    "info_coolant.h5",
    "info_moderator.h5",
    "info_reflector.h5",
    f"{name[0]}.rmc.out",
    f"{name[0]}.rmc.Tally",
    "MeshTally1.h5",
    "MeshTally2.h5",
    f"{name[1]}.dat.h5"
]
destination_folder = os.path.join(script_dir, "histories")
os.makedirs(destination_folder, exist_ok=True)
###########################################################################################

# Step.III, iter1，第一次耦合无法进行松弛
i = 1
os.system(run_MC)
os.system(run_CFD)
for file in source_files:
    file_name, file_ext = os.path.splitext(file)
    new_file_name = f"{file_name}_iter{i}{file_ext}"
    source_path = os.path.join(script_dir, file)
    destination_path = os.path.join(destination_folder, new_file_name)
    shutil.copyfile(source_path, destination_path)
    '''
    if file == "MeshTally1.h5":
        new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
        destination_path = os.path.join(destination_folder, new_file_name)
        shutil.copyfile(source_path, destination_path)
    '''
outfile_path = os.path.join(script_dir, f'{name[0]}.rmc.out')
with open(outfile_path, 'r') as outfile:
    for line in outfile:
        if 'Final Keff:' in line:
            numbers = re.findall(r'[0-9\.]+', line)
            if len(numbers) >= 2:
                keff.append(float(numbers[0]))
                std.append(float(numbers[1]))
                keff_re_diff.append(1)

# Step.IV, 循环运行中子物理及热工水力程序，文件归档，以及松弛迭代过程
i = 2
while i <= max_iter:

    # run RMC for Meshtally1_iter{i}.h5
    os.system(run_MC)

    # 储存keff信息
    outfile_path = os.path.join(script_dir, f'{name[0]}.rmc.out')
    with open(outfile_path, 'r') as outfile:
        for line in outfile:
            if 'Final Keff:' in line:
                numbers = re.findall(r'[0-9\.]+', line)
                if len(numbers) >= 2:
                    keff.append(float(numbers[0]))
                    std.append(float(numbers[1]))
                    keff_re_diff.append(abs(keff[i-1]-keff[i-2])/keff[i-2])

    # 松弛前先进行功率的文件归档
    power_profile = "MeshTally1.h5"
    file_name, file_ext = os.path.splitext(power_profile)
    new_file_name = f"{file_name}_iter{i}{file_ext}"
    source_path = os.path.join(script_dir, power_profile)
    destination_path = os.path.join(destination_folder, new_file_name)
    shutil.copyfile(source_path, destination_path)

    # 由蒙卡结果的标准差计算判敛准则，只读取一个Tally包含的所有行即可
    Tally_path = os.path.join(destination_folder, f'{name[0]}.rmc_iter{i-1}.Tally')
    with open(Tally_path, 'r') as tallyfile:
        line = tallyfile.readline()
        while line:
            if line.startswith('------------------ ID = 1,  Type = power, Number of mesh grids'):
                line = tallyfile.readline()
                line = tallyfile.readline()
                power_re = []
                while not (line.isspace() or line.startswith('====')):
                    power_re.append(float(line.split()[-1]))
                    line = tallyfile.readline()
                break
            line = tallyfile.readline()
        power_re_np = np.array(power_re)

    re_ave.append(np.sum(power_re_np) / np.count_nonzero(power_re_np) * re_factor)

    # 计算方均根，如果没收敛就松弛，否则无需任何处理
    previous_file_name = f"{file_name}_iter{i-1}{file_ext}"
    previous_file_path = os.path.join(destination_folder, previous_file_name)
    with h5py.File(source_path, "r+") as source_file:
        data = source_file["Type2"][:]
        with h5py.File(previous_file_path, "r") as previous_file:
            previous_data = previous_file["Type2"][:]
            non_zero_indices = np.nonzero(previous_data)
            relative_difference = np.divide(data[non_zero_indices], previous_data[non_zero_indices]) - 1
            rms.append(np.sqrt(np.mean(np.square(relative_difference))))
            if ((all_iter_flag == 0) and (rms[i-2] > re_ave[i-2])) or (all_iter_flag == 1):
                updated_data = Omega * data + (1 - Omega) * previous_data
                source_file["Type2"][:] = updated_data

    if ((all_iter_flag == 0) and (rms[i-2] > re_ave[i-2])) or (all_iter_flag == 1):
        new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
        destination_path = os.path.join(destination_folder, new_file_name)
        shutil.copyfile(source_path, destination_path)

    # run fluent for info_fuel_iter{i}.h5 .etc
    os.system(run_CFD)

    # 热工水力计算结果的文件归档与松弛处理判定
    for file in source_files:
        if file == "MeshTally1.h5":
            continue

        file_name, file_ext = os.path.splitext(file)
        new_file_name = f"{file_name}_iter{i}{file_ext}"
        source_path = os.path.join(script_dir, file)
        destination_path = os.path.join(destination_folder, new_file_name)
        shutil.copyfile(source_path, destination_path)

        if (SOR_flag_for_th == 1) and (((all_iter_flag == 0) and (rms[i-2] > re_ave[i-2])) or (all_iter_flag == 1)):
            # 更新 info_fuel.h5
            if file == "info_fuel.h5":
                previous_file_name = f"{file_name}_iter{i-1}{file_ext}"
                previous_file_path = os.path.join(destination_folder, previous_file_name)

                with h5py.File(source_path, "r+") as source_file:
                    data = source_file["temp_fuel"][:]
                    if flag == 1:
                        data1 = source_file["temp_fuel1"][:]
                        data2 = source_file["temp_fuel2"][:]
                        data3 = source_file["temp_fuel3"][:]
                        data4 = source_file["temp_fuel4"][:]
                        data5 = source_file["temp_fuel5"][:]
                    # 进行超松弛迭代更新
                    with h5py.File(previous_file_path, "r") as previous_file:
                        previous_data = previous_file["temp_fuel"][:]
                        updated_data = Omega * data + (1 - Omega) * previous_data
                        if flag == 1:
                            previous_data1 = previous_file["temp_fuel1"][:]
                            previous_data2 = previous_file["temp_fuel2"][:]
                            previous_data3 = previous_file["temp_fuel3"][:]
                            previous_data4 = previous_file["temp_fuel4"][:]
                            previous_data5 = previous_file["temp_fuel5"][:]
                            updated_data1 = Omega * data1 + (1 - Omega) * previous_data1
                            updated_data2 = Omega * data2 + (1 - Omega) * previous_data2
                            updated_data3 = Omega * data3 + (1 - Omega) * previous_data3
                            updated_data4 = Omega * data4 + (1 - Omega) * previous_data4
                            updated_data5 = Omega * data5 + (1 - Omega) * previous_data5
                        # 将更新后的数据写入目标文件
                        source_file["temp_fuel"][:] = updated_data
                        if flag == 1:
                            source_file["temp_fuel1"][:] = updated_data1
                            source_file["temp_fuel2"][:] = updated_data2
                            source_file["temp_fuel3"][:] = updated_data3
                            source_file["temp_fuel4"][:] = updated_data4
                            source_file["temp_fuel5"][:] = updated_data5

                new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
                destination_path = os.path.join(destination_folder, new_file_name)
                shutil.copyfile(source_path, destination_path)

            # 更新 info_coolant.h5
            if file == "info_coolant.h5":
                previous_file_name = f"{file_name}_iter{i-1}{file_ext}"
                previous_file_path = os.path.join(destination_folder, previous_file_name)

                with h5py.File(source_path, "r+") as source_file:
                    data1 = source_file["temp_coolant"][:]
                    data2 = source_file["r_coolant"][:]
                    # 进行超松弛迭代更新
                    with h5py.File(previous_file_path, "r") as previous_file:
                        previous_data1 = previous_file["temp_coolant"][:]
                        previous_data2 = previous_file["r_coolant"][:]
                        updated_data1 = Omega * data1 + (1 - Omega) * previous_data1
                        updated_data2 = Omega * data2 + (1 - Omega) * previous_data2
                        # 将更新后的数据写入目标文件
                        source_file["temp_coolant"][:] = updated_data1
                        source_file["r_coolant"][:] = updated_data2

                new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
                destination_path = os.path.join(destination_folder, new_file_name)
                shutil.copyfile(source_path, destination_path)

            # 更新 info_moderator.h5
            if file == "info_moderator.h5":
                previous_file_name = f"{file_name}_iter{i-1}{file_ext}"
                previous_file_path = os.path.join(destination_folder, previous_file_name)

                with h5py.File(source_path, "r+") as source_file:
                    data = source_file["temp_moderator"][:]
                    # 进行超松弛迭代更新
                    with h5py.File(previous_file_path, "r") as previous_file:
                        previous_data = previous_file["temp_moderator"][:]
                        updated_data = Omega * data + (1 - Omega) * previous_data
                        # 将更新后的数据写入目标文件
                        source_file["temp_moderator"][:] = updated_data

                new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
                destination_path = os.path.join(destination_folder, new_file_name)
                shutil.copyfile(source_path, destination_path)

            # 更新 info_reflector.h5
            if file == "info_reflector.h5":
                previous_file_name = f"{file_name}_iter{i-1}{file_ext}"
                previous_file_path = os.path.join(destination_folder, previous_file_name)

                with h5py.File(source_path, "r+") as source_file:
                    data = source_file["temp_reflector"][:]
                    # 进行超松弛迭代更新
                    with h5py.File(previous_file_path, "r") as previous_file:
                        previous_data = previous_file["temp_reflector"][:]
                        updated_data = Omega * data + (1 - Omega) * previous_data
                        # 将更新后的数据写入目标文件
                        source_file["temp_reflector"][:] = updated_data

                new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
                destination_path = os.path.join(destination_folder, new_file_name)
                shutil.copyfile(source_path, destination_path)

    i += 1
    if (all_iter_flag == 0) and (rms[i-2] <= re_ave[i-2]):
        break


# Step.V, 输出残差信息
save_path = os.path.join(destination_folder, 'residual.dat')
iter_list = np.array([x for x in range(2, i)]).astype(int)
data = np.column_stack((iter_list, rms, re_ave))
iter_list = np.array([x for x in range(1, i)]).astype(int)
keff_data = np.column_stack((iter_list, keff, std, keff_re_diff))
with open(save_path, 'w') as f:
    f.write(f'Omega = {Omega}\n')
    f.write(f'\nSOR for th calculation is [{"OFF" if not (SOR_flag_for_th) else "ON"}]\n')
    f.write(f'\nMultilevel calculation is [{"OFF" if not (flag) else "ON"}]\n')
    if all_iter_flag == 1:
        print(f"\nAll {max_iter} iterations calculated forcefully\n")
        f.write(f"\nAll {max_iter} iterations calculated forcefully\n")
    else:
        if i <= max_iter:
            print(f"\nConverged at iteration {i}\n")
            f.write(f"\nConverged at iteration {i}\n")
        else:
            print(f"\nHave not converged until iter {max_iter}\n")
            f.write(f"\nHave not converged until iter {max_iter}\n")
    f.write('\niter        rms(P(i)/P(i-1)-1)        re_average(P(i-1))\n')
    np.savetxt(f, data, fmt='%d              %.4e                 %.4e')
    f.write('\n\niter    keff      std       keff_re_difference\n')
    np.savetxt(f, keff_data, fmt='%d          %.6f    %.6f    %.6f')

source_path = os.path.join(script_dir, 'coupling.dat')
destination_path = os.path.join(destination_folder, 'coupling.dat')
shutil.copyfile(source_path, destination_path)

os.system('rm -f *.trn')
os.system('ls *.h5 | grep -v ".cas.h5$" | xargs rm -f')
os.system(f'ls {name[0]}* | egrep -v "(.rmc$|.cas.h5$)" | xargs rm -f')
