import os
import shutil
import h5py
import numpy as np

# Step.0, 获取当前脚本所在位置
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)


# Step.I, 读取 coupling.dat 中的迭代信息
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
    max_iter = int(lines[0])


# Step.II, 初始化命令与残差向量，并设置源文件和目标路径
run_MC = f'mpiexec -n {processes_RMC} RMC Singlepin.rmc'
run_CFD = f'fluent 3ddp -g -i runFluent.jou -t{processes_CFD}'
rms = []
rms_SOR = []
re = []

source_files = [
    "info_fuel.h5",
    "info_coolant.h5",
    "info_moderator.h5",
    "info_reflector.h5",
    "Singlepin.rmc.out",
    "Singlepin.rmc.Tally",
    "MeshTally1.h5",
    "MeshTally2.h5",
    "pin.dat.h5"
]
destination_folder = os.path.join(script_dir, "histories")
os.makedirs(destination_folder, exist_ok=True)


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
    if file == "MeshTally1.h5":
        new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
        destination_path = os.path.join(destination_folder, new_file_name)
        shutil.copyfile(source_path, destination_path)


# Step.IV, 循环运行中子物理及热工水力程序，文件归档，以及松弛迭代过程
i = 2
while i <= max_iter:

    # run RMC for Meshtally1_iter{i}.h5
    os.system(run_MC)

    # 松弛前先进行功率的文件归档
    power_profile = "MeshTally1.h5"
    file_name, file_ext = os.path.splitext(power_profile)
    new_file_name = f"{file_name}_iter{i}{file_ext}"
    source_path = os.path.join(script_dir, power_profile)
    destination_path = os.path.join(destination_folder, new_file_name)
    shutil.copyfile(source_path, destination_path)

    # 由蒙卡结果的标准差计算判敛准则，只读取一个Tally包含的所有行即可
    Tally_path = os.path.join(destination_folder, f'Singlepin.rmc_iter{i-1}.Tally')
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
        re_average = np.sum(power_re_np) / np.count_nonzero(power_re_np)
    re.append(re_average * re_factor)

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
            if ((all_iter_flag == 0) and (rms[i-2] > re[i-2])) or (all_iter_flag == 1):
                updated_data = Omega * data + (1 - Omega) * previous_data
                source_file["Type2"][:] = updated_data

    if ((all_iter_flag == 0) and (rms[i-2] > re[i-2])) or (all_iter_flag == 1):
        new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
        destination_path = os.path.join(destination_folder, new_file_name)
        shutil.copyfile(source_path, destination_path)

    # 这段选择性保留！！！这一次仅仅为了观察p与p*哪个量收敛的快，取决于结果可能需要大幅修改源码！如果Omega = 1，则rms向量与rms_SOR向量一致！[这段大概率不需要]
    previous_file_name = f"{file_name}_iter{i-1}_SOR{file_ext}"
    previous_file_path = os.path.join(destination_folder, previous_file_name)
    with h5py.File(destination_path, "r") as source_file:
        data_SOR = source_file["Type2"][:]
        with h5py.File(previous_file_path, "r") as previous_file:
            previous_data_SOR = previous_file["Type2"][:]
            non_zero_indices = np.nonzero(previous_data_SOR)
            relative_difference = np.divide(data_SOR[non_zero_indices], previous_data_SOR[non_zero_indices]) - 1
            rms_SOR.append(np.sqrt(np.mean(np.square(relative_difference))))

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

        if (SOR_flag_for_th == 1) and (((all_iter_flag == 0) and (rms[i-2] > re[i-2])) or (all_iter_flag == 1)):
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

    if (all_iter_flag == 0) and (rms[i-2] <= re[i-2]):
        break

    # 增加计数器
    i += 1


# Step.V, 输出残差信息
save_path = os.path.join(destination_folder, 'residual.dat')
data = np.vstack((rms, rms_SOR, re)).T
header = 'rms          rms_SOR          re'
with open(save_path, 'w') as f:
    f.write(f'Omega = {Omega},   th_SOR = {SOR_flag_for_th}')
    f.write(header + '\n')
    np.savetxt(f, data)

source_path = os.path.join(script_dir, 'coupling.dat')
destination_path = os.path.join(destination_folder, 'coupling.dat')
shutil.copyfile(source_path, destination_path)
