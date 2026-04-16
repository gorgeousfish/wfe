// recompile_mlib.do — Recompile lwfe.mlib from current mata sources
// Can be run from anywhere: stata -b do /path/to/recompile_mlib.do

clear all
set matastrict on

cd "/Users/cxy/Desktop/2026project/sjwfe/wfe-main"

// Source all mata files in dependency order
quietly do mata/wfe_utils.mata
quietly do mata/wfe_ols.mata
quietly do mata/wfe_wwdemean.mata
quietly do mata/wfe_weights_unit.mata
quietly do mata/wfe_weights_time.mata
quietly do mata/wfe_weights_fd.mata
quietly do mata/wfe_weights_did.mata
quietly do mata/wfe_weights_mdid.mata
quietly do mata/wfe_transform.mata
quietly do mata/wfe_se_hac.mata
quietly do mata/wfe_se_fe_twoway.mata
quietly do mata/wfe_se_gmm.mata
quietly do mata/wfe_se_pwfe.mata
quietly do mata/wfe_white_test.mata
quietly do mata/wfe_complex_project.mata
quietly do mata/wfe_twoway_fe_ols.mata
quietly do mata/wfe_postestimation.mata
quietly do mata/wfe_bridge.mata
quietly do mata/wfe_oneway.mata
quietly do mata/wfe_twoway.mata
quietly do mata/wfe_pwfe.mata

// Patch the pscore bridge (overlong identifier)
tempfile __pwfe_bridge_patched
quietly filefilter "mata/wfe_pscore_bridge.mata" "`__pwfe_bridge_patched'", ///
    from("_wfe_pwfe_bridge_numeric_singlevar") ///
    to("_wfe_pwfe_bridge_numvar") replace
quietly do "`__pwfe_bridge_patched'"

// Create the mlib
capture erase lwfe.mlib
mata: mata mlib create lwfe, dir(".") replace
mata: mata mlib add lwfe *()

display as text "lwfe.mlib recompiled successfully"
