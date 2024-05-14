import os
import shutil
import h5py
import numpy as np

# 获取当前脚本所在位置
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# 初始化命令与残差向量
run_MC = 'mpiexec -n 100 RMC Singlepin.rmc'
run_CFD = 'fluent 3ddp -g -i runFluent.jou -t50'
rms = []

# 读取 coupling.dat 中的网格信息
coupling_path = os.path.join(script_dir, 'coupling.dat')
with open(coupling_path, 'r') as coupling_file:
    lines = coupling_file.readlines()
    for i, line in enumerate(lines):
        if line.strip() == '#DIM':
            dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMF':
            out_F_dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMC':
            out_C_dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMM':
            out_M_dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMR':
            out_R_dim = lines[i+1].split()
        elif line.strip() == '#P_bndry in [cm] {xmin xmax ymin ymax zmin zmax}':
            P_bndry = lines[i+1].split()
        elif line.strip() == '#F_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            F_bndry = lines[i+1].split()
        elif line.strip() == '#C_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            C_bndry = lines[i+1].split()
        elif line.strip() == '#M_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            M_bndry = lines[i+1].split()
        elif line.strip() == '#R_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            R_bndry = lines[i+1].split()
        elif line.startswith('#Multilevel Flag'):
            flag = int(lines[i].split()[-1])
        elif line.startswith('#Omega'):
            Omega = float(lines[i].split()[-1])
        elif line.startswith('#Residual'):
            re = float(lines[i].split()[-1])
    max_iter = int(lines[0])

# 设置源文件和目标路径
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

# iter1，第一次耦合无法进行松弛
i = 1

os.system(run_MC)  # run RMC for Meshtally1_iter1.h5
    
os.system(run_CFD)  # run Fluent for info_fuel_iter1.h5 .etc

for file in source_files:    
    file_name, file_ext = os.path.splitext(file)
    new_file_name = f"{file_name}_iter{i}{file_ext}"
    source_path = os.path.join(script_dir, file)
    destination_path = os.path.join(destination_folder, new_file_name)
    shutil.copyfile(source_path, destination_path)

# 循环复制并重命名所有文件，以及松弛迭代过程
i = 2
while i <= max_iter:
    
    os.system(run_MC)  # run RMC for Meshtally1_iter{i}.h5
    
    power_profile = "MeshTally1.h5"
    file_name, file_ext = os.path.splitext(power_profile)
    new_file_name = f"{file_name}_iter{i}{file_ext}"
    source_path = os.path.join(script_dir, power_profile)
    destination_path = os.path.join(destination_folder, new_file_name)
    shutil.copyfile(source_path, destination_path)
    
    # 计算残差，如果没收敛就松弛，否则无需任何处理，松弛前先进行功率的文件归档
    previous_file_name = f"{file_name}_iter{i-1}{file_ext}"
    previous_file_path = os.path.join(destination_folder, previous_file_name)
    with h5py.File(source_path, "r+") as source_file:
        data = source_file["Type2"][:]
        with h5py.File(previous_file_path, "r") as previous_file:
            previous_data = previous_file["Type2"][:]
            non_zero_indices = np.nonzero(previous_data)
            relative_difference = np.divide(data[non_zero_indices], previous_data[non_zero_indices]) - 1
            rms.append(np.sqrt(np.mean(np.square(relative_difference))))
            if rms[i-2] > re:
                updated_data = Omega * data + (1 - Omega) * previous_data
                source_file["Type2"][:] = updated_data

    new_file_name = f"{file_name}_iter{i}_SOR{file_ext}"
    destination_path = os.path.join(destination_folder, new_file_name)
    shutil.copyfile(source_path, destination_path)

    os.system(run_CFD)  # 运行Fluent得到info_fuel_iter{i}.h5
    
    for file in source_files:
        if file == "MeshTally1.h5":
            continue
        
        file_name, file_ext = os.path.splitext(file)
        new_file_name = f"{file_name}_iter{i}{file_ext}"
        source_path = os.path.join(script_dir, file)
        destination_path = os.path.join(destination_folder, new_file_name)
        shutil.copyfile(source_path, destination_path)

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

    if rms[i-2] <= re:
        break
        
    # 增加计数器
    i += 1
    
save_path = os.path.join(destination_folder, 'residual.dat')
np.savetxt(save_path, rms)
