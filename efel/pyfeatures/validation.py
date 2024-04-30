"""Contains scientific validation methods on input signals."""
import efel


def check_ais_initiation(soma_trace: dict, ais_trace: dict) -> bool:
    """Checks the initiation of action potential in AIS with respect to soma."""
    f_values = efel.get_feature_values(
        traces=[soma_trace, ais_trace], feature_names=["AP_begin_time"],
        parallel_map=None, return_list=True, raise_warnings=True
    )
    soma_ap_begin_time = f_values[0]["AP_begin_time"]
    ais_ap_begin_time = f_values[1]["AP_begin_time"]

    if soma_ap_begin_time is None or ais_ap_begin_time is None:
        raise ValueError("AP_begin_time feature is not computed")

    if len(soma_ap_begin_time) != len(ais_ap_begin_time):
        raise ValueError("Number of APs in soma and AIS do not match")

    # Checking if no spike in the soma starts earlier than in the AIS
    for soma_time, ais_time in zip(soma_ap_begin_time, ais_ap_begin_time):
        if soma_time < ais_time:
            raise ValueError(
                "There is a spike that initiates in the soma before the axon"
            )

    return True
