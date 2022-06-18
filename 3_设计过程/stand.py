stand = {
    "K_c": {
        "I":{ # 组合I 无水位 永久+主可变
            "gamma_E1":1.25,
            "gamma_G":0.75,
            "gamma_E2":0.85
        },
        "IV":{ # 组合IV I+地震
            "gamma_E1":1.3,
            "gamma_E2":0.9,
            "gamma_G":0.85
        }
    },
    "K_0":{
        "I":{ # 组合I 无水位 永久+主可变
            "gamma_E1":1.5,
            "gamma_G":1,
            "gamma_E2":0.85
        },
        "IV":{ # 组合IV I+地震
            "gamma_E1":1.35,
            "gamma_E2":1,
            "gamma_G":0.9
        }
    },
    "e":{
        "a":{
            "normal":0.25,
            "earthquake":0.3333333
        },
        "b":{
            "normal":0.16666666,
            "earthquake":0.25
        },
        "c":{
            "normal":0.16666666,
            "earthquake":0.2
        },
        "d":{
            "normal":0.16666666,
            "earthquake":0.16666666
        }
    },
    "sigma":{
        "I":{ # 组合I 永久+主可变
            "gamma_toe":1, #墙趾处地基承载力
            "gamma_heel":1.3, #墙踵处地基承载力
            "gamma_a":1  #地基平均承载力
        },
        "IV":{ # 组合IV I+地震
            "a":{
                "gamma_toe":1.5, #墙趾处地基承载力
                "gamma_heel":1.95, #墙踵处地基承载力
                "gamma_a":1.5  #地基平均承载力
            },
            "b":{
                "gamma_toe":1.4, #墙趾处地基承载力
                "gamma_heel":1.82, #墙踵处地基承载力
                "gamma_a":1.4  #地基平均承载力
            },
            "c":{
                "gamma_toe":1.3, #墙趾处地基承载力
                "gamma_heel":1.69, #墙踵处地基承载力
                "gamma_a":1.3  #地基平均承载力
            },
            "d":{
                "gamma_toe":1.2, #墙趾处地基承载力
                "gamma_heel":1.56, #墙踵处地基承载力
                "gamma_a":1.2  #地基平均承载力
            }
        }
    }
}