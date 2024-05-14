# 设置源文件和目标路径
$sourceFiles = @("E:\Research\2023_ss_Coupling\PWRpin\run\info_fuel.h5", `
"E:\Research\2023_ss_Coupling\PWRpin\run\info_coolant.h5", `
"E:\Research\2023_ss_Coupling\PWRpin\run\info_moderator.h5", `
"E:\Research\2023_ss_Coupling\PWRpin\run\Singlepin.rmc.out", `
"E:\Research\2023_ss_Coupling\PWRpin\run\Singlepin.rmc.Tally", `
"E:\Research\2023_ss_Coupling\PWRpin\run\MeshTally1.h5", `
"E:\Research\2023_ss_Coupling\PWRpin\run\MeshTally2.h5", `
"E:\Research\2023_ss_Coupling\PWRpin\run\pin.dat.h5")
$destinationFolder = "E:\Research\2023_ss_Coupling\PWRpin\run\histories"

.\initialize.exe

$newFileName = "info_fuel_iter0.h5"
$destinationPath = Join-Path -Path $destinationFolder -ChildPath $newFileName
Copy-Item -Path "E:\Research\2023_ss_Coupling\PWRpin\run\info_fuel.h5" -Destination $destinationPath

$newFileName = "info_coolant_iter0.h5"
$destinationPath = Join-Path -Path $destinationFolder -ChildPath $newFileName
Copy-Item -Path "E:\Research\2023_ss_Coupling\PWRpin\run\info_coolant.h5" -Destination $destinationPath

$newFileName = "info_moderator_iter0.h5"
$destinationPath = Join-Path -Path $destinationFolder -ChildPath $newFileName
Copy-Item -Path "E:\Research\2023_ss_Coupling\PWRpin\run\info_moderator.h5" -Destination $destinationPath

# 初始化计数器变量
$i = 1

# 循环复制并重命名文件
while ($i -le 8) {

    .\mpiexec -n 6 RMC.exe Singlepin.rmc

    .\Readh5.exe

    fluent 3ddp -g -i runFluent.jou -t4

    .\Writeh5.exe

    foreach ($file in $sourceFiles) {
        # 生成新文件名
        $newFileName = "{0}_iter{1}{2}" -f [System.IO.Path]::GetFileNameWithoutExtension($file), $i, [System.IO.Path]::GetExtension($file)

        # 拼接目标文件路径
        $destinationPath = Join-Path -Path $destinationFolder -ChildPath $newFileName

        # 复制并重命名文件
        Copy-Item -Path $file -Destination $destinationPath
    }

    # 增加计数器
    $i++
}
