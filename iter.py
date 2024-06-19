import os
import shutil
import h5py
import numpy as np
import re

# Step.0. 获取当前脚本所在位置
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)


# Step.I. 读取 coupling.dat 中的迭代信息
coupling_path = os.path.join(script_dir, 'coupling.dat')
with open(coupling_path, 'r') as coupling_file:
    lines = coupling_file.readlines()
    for i, line in enumerate(lines):
        if line.strip() == '#DIM':
            dim = lines[i+1].split()
        elif line.startswith('#Multilevel Flag'):
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
        elif line.startswith('#tmax discrepancy'):
            epsilon_t = float(lines[i].split()[-1])
        elif line.startswith('#Output all iter flag'):
            all_iter_flag = int(lines[i].split()[-1])
        elif line.startswith('#Project Name'):
            name = lines[i+1].split()
        elif line.startswith('#Coupling mode'):
            mode = int(lines[i+1])
    max_iter = int(lines[0])


# Step.II. 初始化命令，并设置目标路径
destination_folder = os.path.join(script_dir, "histories")
os.makedirs(destination_folder, exist_ok=True)
run_MC = f'mpiexec -n {processes_RMC} RMC {name[0]}.rmc'
run_CFD = f'fluent 3ddp -g -i runFluent{mode}.jou -t{processes_CFD}'


# 稳态情形
def steady_state_cal():
    # Step.III. 初始化时间步信息，残差向量并设置源文件信息
    with open(coupling_path, 'r') as coupling_file:
        lines = coupling_file.readlines()
        for i, line in enumerate(lines):
            if line.startswith('#time step'):
                max_ss_inner_iter_step = int(lines[i+1].split()[4])

    jou_path = os.path.join(script_dir, f'runFluent{mode}.jou')  # runFluent0.jou 文件的路径
    temp_file_path = jou_path + '.tmp'  # 创建一个临时文件用于修改
    with open(jou_path, 'r', encoding='utf-8') as jou_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
        for line in jou_file:
            if line.startswith('rc'):
                temp_file.write(f'rc {name[1]}.cas.h5\n')
            elif line.startswith('wd'):
                temp_file.write(f'wd {name[1]}.dat.h5\n')
            elif line.startswith('solve/iterate'):
                temp_file.write(f'solve/iterate {max_ss_inner_iter_step}\n')
            else:
                temp_file.write(line)
    shutil.copyfile(temp_file_path, jou_path)  # 强制替换原始文件
    os.remove(temp_file_path)  # 删除临时文件

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! change this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    end_files = [
        # "info_fuel.h5",
        # "info_coolant.h5",
        # "info_moderator.h5",
        # "info_reflector.h5",
        f"{name[0]}.rmc.out",
        f"{name[0]}.rmc.Tally",
        f"{name[0]}.rmc.State.h5",
        f"{name[0]}.rmc.innerproduct",
        # "MeshTally1.h5",
        "MeshTally2.h5",
        f"{name[1]}.dat.h5"
    ]
    thermal_files = [
        "info_fuel.h5",
        "info_coolant.h5",
        "info_moderator.h5",
        "info_reflector.h5"
    ]
    ###########################################################################################

    i = 1
    converge = False
    n_P_converge = False
    n_k_converge = False
    th_converge = False
    L2_norm = []                    # 【MC分布量】单时间步内两次Picard迭代间三维功率差向量的二范数
    Linf = []                       # 【MC分布量】单时间步内两次Picard迭代间三维功率差向量的无穷范数，比二范数更加严格
    re_ave = []                     # 【MC分布量】前一次Picard功率向量的元素平均标准差
    keff = []                       # 【MC集总量】keff列表
    std = []                        # 【MC集总量】keff的标准差列表
    keff_re_diff = []               # 【MC集总量】keff的相对偏差
    tmax = []                       # 【CFD集总量】最高燃料温度列表，由于标准差难以确定，故不采用分布量作为热工水力判敛条件
    tmax_re_diff = []               # 【CFD集总量】最高燃料温度相对偏差列表

    # Step.IV. 循环运行中子物理及热工水力程序，文件归档，以及松弛迭代过程，第一次耦合无法进行松弛
    while i <= max_iter and (all_iter_flag == 1 or not converge):

        # Step.IV.1. run RMC with post process
        os.system(run_MC)
        n_k_converge, n_P_converge = MC_post(i, destination_folder, keff, std, keff_re_diff, L2_norm, Linf, re_ave)

        # Step.IV.2. run Fluent with post process
        os.system(run_CFD)
        th_converge = CFD_post(i, destination_folder, thermal_files, tmax, tmax_re_diff)

        # Step.IV.3. 当且仅当蒙卡计算的集总量与分布量，CFD计算的集总量全部收敛时，判定耦合收敛
        converge = n_k_converge * n_P_converge * th_converge

        # Step.IV.4. 将其余结果文件归档
        for file in end_files:
            if os.path.exists(os.path.join(script_dir, file)):
                first_dot_index = file.find('.')
                file_name, file_ext = (file[:first_dot_index], file[first_dot_index:]) if first_dot_index != -1 else (file, '')
                shutil.copyfile(os.path.join(script_dir, file), os.path.join(destination_folder, f"{file_name}_iter{i}{file_ext}"))
            else:
                print(f"\nWarning: File {file} does not exist and will be skipped, please check the name.\n")

        i += 1

    # Step.V. 输出残差信息等耦合结果
    Picard_output(i, destination_folder, L2_norm, Linf, re_ave, keff, std, keff_re_diff, tmax, tmax_re_diff)

    # Step.VI. 保存耦合控制输入卡，并清理工作路径的临时文件
    shutil.copyfile(os.path.join(script_dir, 'coupling.dat'), os.path.join(destination_folder, 'coupling.dat'))
    os.system('rm -f *.trn')
    os.system('ls *.h5 | grep -v ".cas.h5$" | xargs rm -f')
    os.system(f'ls {name[0]}* | egrep -v "(.rmc$|.cas.h5$)" | xargs rm -f')


# 瞬态情形
def transient():
    # Step.III. 初始化时间步信息，读取初始功率与功率因子，并设置源文件信息
    with open(coupling_path, 'r') as coupling_file:
        lines = coupling_file.readlines()
        for i, line in enumerate(lines):
            if line.startswith('#time step'):
                deltat = float(lines[i+1].split()[0])
                max_time = float(lines[i+1].split()[1])
                Nneu = int(lines[i+1].split()[2])
                Nth = int(lines[i+1].split()[3])
                max_inner_th_step = int(lines[i+1].split()[4])
                sub_deltat_th = float(deltat/Nth)
            elif line.startswith('#Coupling order'):
                if int(lines[i+1]) == 1:
                    order_flag = 1
                else:
                    order_flag = 0
        init_power = float(lines[1])
    with open(os.path.join(script_dir, f'{name[0]}.rmc.innerproduct'), 'r') as init_power_file:
        lines = init_power_file.readlines()
        for line in lines:
            if line.startswith('Average total power:'):
                init_power_factor = float(line.split(':')[-1])

    rmc_input_path = os.path.join(script_dir, f'{name[0]}.rmc')  # .rmc 文件的路径
    temp_file_path = rmc_input_path + '.tmp'  # 创建一个临时文件用于修改
    with open(rmc_input_path, 'r', encoding='utf-8') as rmc_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
        for line in rmc_file:
            if line.startswith('Timestep deltat'):
                temp_file.write(f'Timestep deltat = {deltat}\n')
            elif line.startswith('Substep number'):
                temp_file.write(f'Substep number = {Nneu} interpolate = 1\n')
            else:
                temp_file.write(line)
    shutil.copyfile(temp_file_path, rmc_input_path)  # 强制替换原始文件
    os.remove(temp_file_path)  # 删除临时文件

    jou_path = os.path.join(script_dir, f'runFluent{mode}.jou')  # runFluent1.jou 文件的路径
    temp_file_path = jou_path + '.tmp'  # 创建一个临时文件用于修改
    with open(jou_path, 'r', encoding='utf-8') as jou_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
        for line in jou_file:
            if line.startswith('rc'):
                temp_file.write(f'rc {name[1]}.cas.h5\n')
            elif line.startswith('rd'):
                temp_file.write(f'rd {name[1]}_init.dat.h5\n')
            elif line.startswith('solve/time-step-size'):
                temp_file.write(f'solve/time-step-size {sub_deltat_th}\n')
            elif line.startswith('solve/dual-time-iterate'):
                temp_file.write(f'solve/dual-time-iterate {Nth} {max_inner_th_step}\n')
            elif line.startswith('wd'):
                temp_file.write(f'wd {name[1]}.dat.h5\n')
            else:
                temp_file.write(line)
    shutil.copyfile(temp_file_path, jou_path)  # 强制替换原始文件
    os.remove(temp_file_path)  # 删除临时文件

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! change this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    init_files = [
        "info_fuel.h5",
        "info_coolant.h5",
        "info_moderator.h5",
        "info_reflector.h5",
        f"{name[0]}_init.rmc.State.h5",  # 这个文件在每个时间步计算结束后需要被该时间步结尾的文件替代
        f"{name[0]}.rmc.innerproduct",
        "MeshTally1.h5",
        f"{name[1]}_init.dat.h5"         # 这个文件在每个时间步计算结束后需要被该时间步结尾的文件替代
    ]

    end_files = [
        # "info_fuel.h5",
        # "info_coolant.h5",
        # "info_moderator.h5",
        # "info_reflector.h5",
        f"{name[0]}.rmc.out",
        f"{name[0]}.rmc.Tally",
        f"{name[0]}.rmc.State.h5",
        f"{name[0]}.rmc.innerproduct",
        # "MeshTally1.h5",
        "MeshTally2.h5",
        f"{name[1]}.dat.h5"
    ]

    thermal_files = [
        "info_fuel.h5",
        "info_coolant.h5",
        "info_moderator.h5",
        "info_reflector.h5"
    ]
    ###########################################################################################

    # Step.IV. 含时的外迭代，循环体外定义的列表储存每个时间步结尾的值，循环体内定义的列表储存每个单时间步内的迭代计算值
    t = 0.0
    n = 1
    time = []
    iteration_step = []
    power_time = []
    power_re_factor_time = []
    keff_time = []
    reactivity_time = []
    tmax_time = []

    time.append(t)
    power_time.append(init_power)
    power_re_factor_time.append(1.0)
    keff_time.append(1.0)
    reactivity_time.append(0.0)
    with h5py.File(os.path.join(script_dir, "info_fuel.h5"), "r") as source_file:
        tmax_time.append(source_file["temp_fuel_max"][()])

    while (t + deltat) <= max_time:
        i = 1  # 重置Picard迭代步计数
        converge = False
        n_P_converge = False
        n_k_converge = False
        th_converge = False
        current_timestep_path = os.path.join(destination_folder, f'timestep_{n}')
        os.makedirs(current_timestep_path, exist_ok=True)
        for file in init_files:
            first_dot_index = file.find('.')
            file_name, file_ext = (file[:first_dot_index], file[first_dot_index:]) if first_dot_index != -1 else (file, '')
            shutil.copyfile(os.path.join(script_dir, file), os.path.join(current_timestep_path, f"{file_name}_iter0{file_ext}"))
        L2_norm = []                    # 【MC分布量】单时间步内两次Picard迭代间三维功率差向量的二范数/root(N) (差向量元素方均根)
        Linf = []                       # 【MC分布量】单时间步内两次Picard迭代间三维功率差向量的无穷范数，比二范数更加严格
        re_ave = []                     # 【MC分布量】前一次Picard功率向量的元素平均标准差
        keff = []                       # 【MC集总量】keff列表
        std = []                        # 【MC集总量】keff的标准差列表
        keff_re_diff = []               # 【MC集总量】keff的相对偏差
        tmax = []                       # 【CFD集总量】最高燃料温度列表，由于标准差难以确定，故不采用分布量作为热工水力判敛条件
        tmax_re_diff = []               # 【CFD集总量】最高燃料温度相对偏差列表

        # Step.V. 单时间步内的Picard内迭代，当max_iter=1时，退化为单步法
        while i <= max_iter and (all_iter_flag == 1 or not converge):

            if order_flag == 0:
                # Step.V.1. run RMC with post process
                shutil.copyfile(os.path.join(script_dir, f'{name[0]}_init.rmc.State.h5'), os.path.join(script_dir, f'{name[0]}.rmc.State.h5'))
                os.system(run_MC)
                n_k_converge, n_P_converge = MC_post(i, current_timestep_path, keff, std, keff_re_diff, L2_norm, Linf, re_ave)  # i与path为输入参数，其余列表作更新
                current_power = transient_power_update(init_power, init_power_factor)  # 计算，储存功率值，并更新coupling.dat以备UDF读取

                # Step.V.2. run Fluent with post process
                os.system(run_CFD)
                th_converge = CFD_post(i, current_timestep_path, thermal_files, tmax, tmax_re_diff)  # i, path与文件列表为输入参数，其余列表作更新

            elif order_flag == 1:
                # Step.V.1. run Fluent with post process
                os.system(run_CFD)
                th_converge = CFD_post(i, current_timestep_path, thermal_files, tmax, tmax_re_diff)  # i, path与文件列表为输入参数，其余列表作更新

                # Step.V.2. run RMC with post process
                shutil.copyfile(os.path.join(script_dir, f'{name[0]}_init.rmc.State.h5'), os.path.join(script_dir, f'{name[0]}.rmc.State.h5'))
                os.system(run_MC)
                n_k_converge, n_P_converge = MC_post(i, current_timestep_path, keff, std, keff_re_diff, L2_norm, Linf, re_ave)  # i与path为输入参数，其余列表作更新
                current_power = transient_power_update(init_power, init_power_factor)  # 计算，储存功率值，并更新coupling.dat以备UDF读取

            # Step.V.3. 当且仅当蒙卡计算的集总量与分布量，CFD计算的集总量全部收敛时，判定耦合收敛
            converge = n_k_converge * n_P_converge * th_converge

            # Step.V.4. 将其余结果文件归档
            for file in end_files:
                if os.path.exists(os.path.join(script_dir, file)):
                    first_dot_index = file.find('.')
                    file_name, file_ext = (file[:first_dot_index], file[first_dot_index:]) if first_dot_index != -1 else (file, '')
                    shutil.copyfile(os.path.join(script_dir, file), os.path.join(current_timestep_path, f"{file_name}_iter{i}{file_ext}"))
                else:
                    print(f"\nWarning: File {file} does not exist and will be skipped, please check the name.\n")

            i += 1

        # Step.VII. 输出单时间步的残差信息等耦合结果，并记录单时间步内最终迭代步的结果以备输出
        Picard_output(i, current_timestep_path, L2_norm, Linf, re_ave, keff, std, keff_re_diff, tmax, tmax_re_diff)

        t += deltat
        n += 1
        time.append(t)
        iteration_step.append(i-1)
        power_time.append(current_power)
        power_re_factor_time.append(current_power/init_power)
        keff_time.append(keff[-1])
        with h5py.File(os.path.join(script_dir, f"{name[0]}.rmc.State.h5"), "r") as source_file:
            reactivity_time.append(source_file["kinetics/reactivity"][()])
        tmax_time.append(tmax[-1])

        # Step.VIII. 将下一个时间步的初始条件替换为本时间步末尾条件，由于后续可能要考虑接续瞬态计算，因此最后一个时间步之后也不宜删除_init文件
        shutil.copyfile(os.path.join(script_dir, f"{name[0]}.rmc.State.h5"), os.path.join(script_dir, f"{name[0]}_init.rmc.State.h5"))
        shutil.copyfile(os.path.join(script_dir, f"{name[1]}.dat.h5"), os.path.join(script_dir, f"{name[1]}_init.dat.h5"))

    # Step.IX. 输出时间相关结果 (to do)
    transient_output(n, destination_folder, time, iteration_step, power_time, power_re_factor_time, keff_time,
                     reactivity_time, tmax_time, deltat, max_time, Nneu, Nth, max_inner_th_step)
    # Step.X. 保存耦合控制输入卡，并清理工作路径的临时文件
    shutil.copyfile(os.path.join(script_dir, 'coupling.dat'), os.path.join(destination_folder, 'coupling.dat'))
    os.system('rm -f *.trn')
    os.system('ls *.h5 | grep -v ".cas.h5$" | xargs rm -f')
    os.system(f'ls {name[0]}* | egrep -v "(.rmc$|.cas.h5$)" | xargs rm -f')


def MC_post(i, path, keff, std, keff_re_diff, L2_norm, Linf, re_ave):
    # 储存keff信息
    outfile_path = os.path.join(script_dir, f'{name[0]}.rmc.out')
    with open(outfile_path, 'r') as outfile:
        for line in outfile:
            if 'Final Keff:' in line:
                numbers = re.findall(r'[0-9\.]+', line)
                if len(numbers) >= 2:
                    keff.append(float(numbers[0]))
                    std.append(float(numbers[1]))
                    if i == 1:
                        keff_re_diff.append(1)
                    else:
                        keff_re_diff.append(abs(keff[i-1]-keff[i-2])/keff[i-2])

    # !!!!!!!!!!!!!!!!!!!!!!  判敛准则  重要  !!!!!!!!!!!!!!!!!!!!!! #
    n_k_converge = keff_re_diff[i-1] < re_factor * std[i-1]
    # !!!!!!!!!!!!!!!!!!!!!!  判敛准则  重要  !!!!!!!!!!!!!!!!!!!!!! #

    # 松弛前先进行功率的文件归档，否则功率原始数据将直接被松弛处理后的数据覆盖，从而损失原始数据
    power_profile = "MeshTally1.h5"
    first_dot_index = power_profile.find('.')
    file_name, file_ext = (power_profile[:first_dot_index], power_profile[first_dot_index:]) if first_dot_index != -1 else (power_profile, '')
    shutil.copyfile(os.path.join(script_dir, power_profile), os.path.join(path, f"{file_name}_iter{i}{file_ext}"))

    if i == 1:
        n_P_converge = False
    elif i > 1:
        # 由蒙卡结果的标准差计算判敛准则，只读取一个Tally包含的所有行即可，注意这句话的含义表明了一定不能是第一次迭代，否则无需也无法计算判敛条件
        Tally_path = os.path.join(path, f'{name[0]}_iter{i-1}.rmc.Tally')
        with open(Tally_path, 'r') as tallyfile:
            line = tallyfile.readline()
            while line:
                if line.startswith('------------------ ID = 1,  Type = power, Number of mesh grids'):
                    line = tallyfile.readline()
                    power_re = []
                    for j in range(0, (int)(dim[0])*(int)(dim[1])*(int)(dim[2])):
                        line = tallyfile.readline()
                        if float(line.split()[-1]) != 0:
                            power_re.append(float(line.split()[-1]))
                    break
                line = tallyfile.readline()
            power_re_np = np.array(power_re)
        re_ave.append(np.sum(power_re_np) / np.count_nonzero(power_re_np) * re_factor)

        # 计算方均根，如果没收敛就[在当前路径]松弛以便CFD求解器直接读取，注意这句话的含义表明了一定不能是第一次迭代，否则无需也无法计算判敛条件
        with h5py.File(os.path.join(script_dir, power_profile), "r+") as source_file:
            data = source_file["Type2"][:]
            with h5py.File(os.path.join(path, f"{file_name}_iter{i-1}{file_ext}"), "r") as previous_file:
                previous_data = previous_file["Type2"][:]
                non_zero_indices = np.nonzero(previous_data)
                abs_difference = np.subtract(data[non_zero_indices], previous_data[non_zero_indices])
                L2_norm.append(np.sqrt(np.mean(np.square(abs_difference))))
                Linf.append(np.max(np.abs(abs_difference)))

                # !!!!!!!!!!!!!!!!!!!!!!  判敛准则  重要  !!!!!!!!!!!!!!!!!!!!!! #
                n_P_converge = max(L2_norm[i-2], Linf[i-2]) <= re_ave[i-2]
                # !!!!!!!!!!!!!!!!!!!!!!  判敛准则  重要  !!!!!!!!!!!!!!!!!!!!!! #

                if all_iter_flag == 1 or not n_P_converge:
                    updated_data = Omega * data + (1 - Omega) * previous_data
                    source_file["Type2"][:] = updated_data

        # 松弛处理后的文件重新归档
        if all_iter_flag == 1 or not n_P_converge:
            shutil.copyfile(os.path.join(script_dir, power_profile), os.path.join(path, f"{file_name}_iter{i}_SOR{file_ext}"))

    return n_k_converge, n_P_converge


def transient_power_update(init_power, init_power_factor):
    # 由.innerproduct提取功率因子的变化并改写耦合输入卡以便UDF读取
    with open(os.path.join(script_dir, f'{name[0]}.rmc.innerproduct'), 'r') as current_power_file:
        lines = current_power_file.readlines()
        for line in lines:
            if line.startswith('Average total power:'):
                current_power_factor = float(line.split(':')[-1])
    current_power = init_power / init_power_factor * current_power_factor
    with open(coupling_path, 'r+') as coupling_file:
        lines = coupling_file.readlines()
        lines[1] = f"{current_power}\n"
        coupling_file.seek(0)
        coupling_file.truncate()
        coupling_file.writelines(lines)
    return current_power


def CFD_post(i, path, thermal_files, tmax, tmax_re_diff):
    # 热工水力计算结果的文件归档与松弛处理判定，同样需要先归档后松弛，以防原始数据被松弛后的结果覆盖，但根据实测，不推荐做松弛
    # 先归档
    for file in thermal_files:
        first_dot_index = file.find('.')
        file_name, file_ext = (file[:first_dot_index], file[first_dot_index:]) if first_dot_index != -1 else (file, '')
        shutil.copyfile(os.path.join(script_dir, file), os.path.join(path, f"{file_name}_iter{i}{file_ext}"))

    # 储存tmax值，并计算判敛条件，注意第一次迭代虽然可以计算，但结果相当于时间步结尾与初始值相对比，无意义，也无需松弛
    with h5py.File(os.path.join(script_dir, "info_fuel.h5"), "r") as source_file:
        current_t = source_file["temp_fuel_max"][()]
        tmax.append(current_t)
        if i == 1:
            tmax_re_diff.append(1)
        else:
            tmax_re_diff.append(abs(tmax[i-1]-tmax[i-2])/tmax[i-2])

    # !!!!!!!!!!!!!!!!!!!!!!  判敛准则  重要  !!!!!!!!!!!!!!!!!!!!!! #
    th_converge = tmax_re_diff[i-1] <= epsilon_t
    # !!!!!!!!!!!!!!!!!!!!!!  判敛准则  重要  !!!!!!!!!!!!!!!!!!!!!! #

    if i > 1:
        # 执行松弛，由于各个.h5文件的数据结构不同，因此是否使用循环体代码长度相仿，此处不用
        if (SOR_flag_for_th == 1) and (all_iter_flag == 1 or not th_converge):
            # 更新 info_fuel.h5
            with h5py.File(os.path.join(script_dir, "info_fuel.h5"), "r+") as source_file:
                data = source_file["temp_fuel"][:]
                if flag == 1:
                    data1 = source_file["temp_fuel1"][:]
                    data2 = source_file["temp_fuel2"][:]
                    data3 = source_file["temp_fuel3"][:]
                    data4 = source_file["temp_fuel4"][:]
                    data5 = source_file["temp_fuel5"][:]
                # 进行超松弛迭代更新
                with h5py.File(os.path.join(path, f"info_fuel_iter{i-1}.h5"), "r") as previous_file:
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
            # 再次归档
            shutil.copyfile(os.path.join(script_dir, "info_fuel.h5"), os.path.join(path, f"info_fuel_iter{i}_SOR.h5"))

            # 更新 info_coolant.h5
            with h5py.File(os.path.join(script_dir, "info_coolant.h5"), "r+") as source_file:
                data1 = source_file["temp_coolant"][:]
                data2 = source_file["r_coolant"][:]
                # 进行超松弛迭代更新
                with h5py.File(os.path.join(path, f"info_coolant_iter{i-1}.h5"), "r") as previous_file:
                    previous_data1 = previous_file["temp_coolant"][:]
                    previous_data2 = previous_file["r_coolant"][:]
                    updated_data1 = Omega * data1 + (1 - Omega) * previous_data1
                    updated_data2 = Omega * data2 + (1 - Omega) * previous_data2
                    # 将更新后的数据写入目标文件
                    source_file["temp_coolant"][:] = updated_data1
                    source_file["r_coolant"][:] = updated_data2
            # 再次归档
            shutil.copyfile(os.path.join(script_dir, "info_coolant.h5"), os.path.join(path, f"info_coolant_iter{i}_SOR.h5"))

            # 更新 info_moderator.h5
            with h5py.File(os.path.join(script_dir, "info_moderator.h5"), "r+") as source_file:
                data = source_file["temp_moderator"][:]
                # 进行超松弛迭代更新
                with h5py.File(os.path.join(path, f"info_moderator_iter{i-1}.h5"), "r") as previous_file:
                    previous_data = previous_file["temp_moderator"][:]
                    updated_data = Omega * data + (1 - Omega) * previous_data
                    # 将更新后的数据写入目标文件
                    source_file["temp_moderator"][:] = updated_data
            # 再次归档
            shutil.copyfile(os.path.join(script_dir, "info_moderator.h5"), os.path.join(path, f"info_moderator_iter{i}_SOR.h5"))

            # 更新 info_reflector.h5
            with h5py.File(os.path.join(script_dir, "info_reflector.h5"), "r+") as source_file:
                data = source_file["temp_reflector"][:]
                # 进行超松弛迭代更新
                with h5py.File(os.path.join(path, f"info_fuel_iter{i-1}.h5"), "r") as previous_file:
                    previous_data = previous_file["temp_reflector"][:]
                    updated_data = Omega * data + (1 - Omega) * previous_data
                    # 将更新后的数据写入目标文件
                    source_file["temp_reflector"][:] = updated_data
            # 再次归档
            shutil.copyfile(os.path.join(script_dir, "info_reflector.h5"), os.path.join(path, f"info_reflector_iter{i}_SOR.h5"))

    return th_converge


def Picard_output(i, path, L2_norm, Linf, re_ave, keff, std, keff_re_diff, tmax, tmax_re_diff):
    save_path = os.path.join(path, 'results.dat')
    iter_list = np.array([x for x in range(2, i)]).astype(int)
    P_data = np.column_stack((iter_list, L2_norm, Linf, re_ave))
    iter_list = np.array([x for x in range(1, i)]).astype(int)
    keff_data = np.column_stack((iter_list, keff, std, keff_re_diff))
    tmax_data = np.column_stack((iter_list, tmax, tmax_re_diff))
    with open(save_path, 'w') as f:
        f.write(f'\nCalculation mode is [{"Steady state" if not (mode) else "Transient"}]\n')
        f.write(f'\nOmega = {Omega}\n')
        f.write(f'\nSOR for th calculation is [{"OFF" if not (SOR_flag_for_th) else "ON"}]\n')
        f.write(f'\nMultilevel calculation is [{"OFF" if not (flag) else "ON"}]\n')
        if all_iter_flag == 1:
            print(f"\nAll {max_iter} iterations are forced to be calculated\n")
            f.write(f"\nAll {max_iter} iterations are forced to be calculated\n")
        else:
            if i <= max_iter:
                print(f"\nConverged at iteration {i}\n")
                f.write(f"\nConverged at iteration {i}\n")
            else:
                print(f"\nHave not converged until iter {max_iter}\n")
                f.write(f"\nHave not converged until iter {max_iter}\n")
        f.write('\niter        L2[P(i)-P(i-1)]/root(N)        Linf[P(i)-P(i-1)]        re_average(P(i-1))\n')
        np.savetxt(f, P_data, fmt='%d              %.4e                 %.4e                 %.4e')
        f.write('\n\niter    keff      std       keff_re_difference\n')
        np.savetxt(f, keff_data, fmt='%d          %.6f    %.6f    %.6f')
        f.write('\n\niter    tmax      tmax_re_difference\n')
        np.savetxt(f, tmax_data, fmt='%d          %.4f    %.4f')


def transient_output(n, path, time, iteration_step, power_time, power_re_factor_time, keff_time,
                     reactivity_time, tmax_time, deltat, max_time, Nneu, Nth, max_inner_th_step):
    save_path = os.path.join(path, 'transient_results.dat')
    iter_list = np.array([x for x in range(0, n)]).astype(int)
    transient_data = np.column_stack((iter_list, time, iteration_step, power_time, power_re_factor_time, keff_time, reactivity_time, tmax_time))
    with open(save_path, 'w') as f:
        f.write('\nTransient results\n')
        f.write(f'\nTime step size = {deltat} s, maximum time = {max_time} s\n')
        f.write(f'\nMC subtimestep number = {Nneu}, CFD subtimestep number = {Nth}, maximum iteration per CFD subtimestep = {max_inner_th_step}\n')
        f.write('\nstep        time        Converge at Picard iteration        power, MW\
        power factor        keff        reactivity        maximum temperature, K\n')
        np.savetxt(f, transient_data, fmt='%d    %.4f    %d    %.6f    %.6f    %.6f    %.6f    %.4f')


if mode == 0:
    steady_state_cal()
elif mode == 1:
    transient()
else:
    print("mode definition error! mode 0 for steady-state, mode 1 for transient")
