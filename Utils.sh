#!/bin/bash

edit_and_run() {
    local file_name_conf="$1"
    local file_name_init="$2"
    local file_number="$3"
    local number_of_files="$4"
    local number_of_events="$5"

    for ((i=0; i<number_of_files; i++)); do
        # Edit the last line of the file
        sed -i '' "$ s/\([^0-9]*\)[0-9]*$/\1$((file_number+i))/" "$file_name_conf"

        # Run the command
        /Users/fabian/Documents/nexus/build/nexus -n "$number_of_events" "$file_name_init"
    done
}


generate_macros() {
    local lengths=("$@")

    # Directory and base config file
    macro_dir="macros"
    base_single_particle_macro="${macro_dir}/generators/SingleParticle.mac"
    base_init_macro="${macro_dir}/cyl_eg.init.mac"
    base_config_macro="${macro_dir}/cyl_eg.config.mac"

    # Check if the base config file exists
    if [ ! -f $base_config_macro ]; then
        echo "Base config macro '$base_config_macro' not found"
        return 1
    fi

    # Check if the base single particle macro exists
    if [ ! -f $base_single_particle_macro ]; then
        echo "Base single particle macro '$base_single_particle_macro' not found"
        return 1
    fi

    # Check if the base init file exists
    if [ ! -f $base_init_macro ]; then
        echo "Base init macro '$base_init_macro' not found"
        return 1
    fi

    # Loop through fiber lengths and create the corresponding macro files
    for i in "${lengths[@]}"; do
        echo $i
        new_config_macro="macros_temp/cyl_eg_${i}.config.mac"
        new_init_macro="macros_temp/cyl_eg_${i}.init.mac"
        new_single_particle_macro="macros_temp/SingleParticle_${i}.mac"
        
        # Copy the base files
        cp "$base_config_macro" "$new_config_macro"
        cp "$base_single_particle_macro" "$new_single_particle_macro"
        cp "$base_init_macro" "$new_init_macro"

        # Replace the last line with the length in the single particle macro
        sed -i '' "s|/source/geometries/CylindricChamber/set_vertex_z .*|/source/geometries/CylindricChamber/set_vertex_z ${i} mm|" "$new_single_particle_macro"

        # New: Replaces the specific pattern wherever it's found
        sed -i '' "s|/control/execute macros/generators/SingleParticle.mac|/control/execute ${new_single_particle_macro}|" "$new_config_macro"

        # Replace the last line in the config macro with the output filename
        sed -i '' "$ s|.*| \/nexus\/persistency\/output_file outputs\/NewCS\/1mil_events_1350\/${i}.next|" "$new_config_macro"
        # 12 means 1.2 m of attenuation length set in the optical material properties

        # Replace the last line in the init macro with the corresponding config macro
        sed -i '' "$ s|.*| \/nexus\/RegisterMacro ${new_config_macro}|" "$new_init_macro"
        
        echo "Created config file: $new_config_macro"
        echo "Created init file: $new_init_macro"
        echo "Created single particle macro: $new_single_particle_macro"
        
    done
}