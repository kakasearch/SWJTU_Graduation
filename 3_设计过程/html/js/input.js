var submit_data = {
    "design":{},
    "tu":{},
    "parms":{}
}
var scale = 30
var offset =  [80,550] 
var safe =  `<svg t="1650422551299" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2862" width="16" height="16"><path d="M68 528.4s229.6 188 252.4 347.6c0 0 301.2-492.8 635.6-536.8 0 0-102-74.4-68.4-191.6 0 0-185.6 18.4-535.2 561.2l-164-278.4L68 528.4z m0 0" fill="#22AC38" p-id="2863"></path></svg>`
var not_safe = `<svg t="1650423738932" class="icon" viewBox="0 0 1025 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="9073" width="16" height="16"><path d="M718.882684 511.351282 1010.410118 800.927611C1024.685139 815.107052 1027.762671 835.153596 1017.28356 845.70286L849.127929 1014.985224C838.649124 1025.534489 818.582055 1022.591441 804.307034 1008.412L512.781438 718.837509 223.267604 1010.296628C209.089387 1024.570118 189.044682 1027.647038 178.496643 1017.169458L9.230209 849.03282C-1.318137 838.554934 1.624604 818.49001 15.802821 804.216827L305.314511 512.759546 13.70467 223.101423C-0.570351 208.921981-3.647577 188.875438 6.831534 178.326173L174.986859 9.043809C185.46597-1.505455 205.533039 1.437592 219.80806 15.617034L511.415756 305.273319 801.039568 13.703302C815.217785-0.570187 835.26249-3.647107 845.810529 6.83078L1015.076963 174.967417C1025.625309 185.444997 1022.682568 205.509921 1008.504351 219.78341L718.882684 511.351282Z" p-id="9074" fill="#d81e06"></path></svg>`
var check_flag = true
function change_BT(){
    if (document.querySelector("#have_tenon").checked){
        for(let i of ["#HT","#BT","#BT1"]){
            document.querySelector(i).disabled = false
        }
    }else{
        for(let i of ["#HT","#BT","#BT1"]){
            document.querySelector(i).disabled = true
        }
    }
}
function check_limit(e,min_,max_){
    //检查输入是否超出限制
    if (e.value > max_ || e.value < min_){
        alert(`${e.parentElement.innerText} 超出规范限制(${min_},${max_})，请检查数据`)
    }
}

//保留n位小数
function roundFun(value, n) {
    return Math.round(value*Math.pow(10,n))/Math.pow(10,n);
  }
function get_design_data(){
    let Decimal = 5 //小数点后几位
    let design = {
        //用户输入
		"H" : 8.5, //挡土墙高度 
		"H2" : 0.5, //趾板高 
		"H3" : 0.5, //踵板高
		"B" : 0.4, //立壁板顶宽   
		"B2" : 1.5, //趾板宽 
		"B3" : 2.5, //踵板宽
		"Hzd" : 0.75, //趾端高度
		"L" : 6,//每节扶壁墙长度
		"B_fb" : 0.5, //扶壁宽度
		"num_fb":2,//单节扶壁墙扶壁数量
		"have_tenon" : true, //是否设凸榫
		"HT" : 0.5, //凸榫高度 
		"BT" : 0.5, //凸榫宽度 
		"BT1" : 2, //凸榫外缘距墙趾的距离
        "cover_h":1.5,//墙趾埋深
    }
    for(let i of ["H","H2","H3","B","B2","B3","Hzd","L","B_fb","num_fb","cover_h"]){
        let v = document.getElementById(i).value
        if(v){
            design[i] = roundFun(parseFloat(v),Decimal)
        }
    }
    if(document.getElementById("have_tenon").checked){
        design["have_tenon"] = true
        for(let i of ["HT","BT","BT1"]){
            let v = document.getElementById(i).value
            if(v){
                design[i] =  roundFun( parseFloat(v),Decimal)
            }
        }
    }else{
        design["have_tenon"] = false
        design["HT"] = 0
        design["BT"] = 0
        design["BT1"] = 0
    }
    //以墙趾处为原点
    let H = design["H"]
    let H2 = design["H2"]
    let H3 = design["H3"]
    let B = design["B"]
    let B2 = design["B2"]
    let B3 = design["B3"]
    let Hzd = design["Hzd"]
    let HT = design["HT"]
    let BT = design["BT"]
    let BT1 = design["BT1"]

    let coords =  [[0,0],[0,H2],[B2,Hzd],[B2,H],[B2+B,H],[B2+B,H3],[B2+B+B3,H3],[B2+B+B3,0]]
    if(design["have_tenon"]){
        coords.push([BT1+BT,0],[BT1+BT,-HT],[BT1,-HT],[BT1,0],[0,0])
    }else{
        coords.push([0,0])
    }
    for(i of coords){
        i[0] = roundFun(i[0],Decimal)
        i[1] = roundFun(i[1],Decimal)
    }
    design["coords"] = coords
    design["coord_qz"] = [ roundFun(B+B2+B3,Decimal),0]
    design["coord_qd"] = [ roundFun(B2+B,Decimal),H]
    submit_data["design"] = design
}
function clearCanvas(canvas_id)
{  
    var c=document.getElementById(canvas_id);
    c.height=c.height;  
}  

function polygon(canvas_id,coords) {
    //绘制多边形
    var canvas = document.getElementById(canvas_id);
    var context = canvas.getContext("2d");
    context.strokeStyle = "#008"; 
    context.lineWidth = 2; // 2个像素宽
    context.beginPath();
    context.moveTo(offset[0]+coords[0][0]*scale,offset[1] -coords[0][1]*scale);
    for (let i=1; i<coords.length; i++) {
        context.lineTo(offset[0]+coords[i][0]*scale,offset[1]-coords[i][1]*scale);
    }
    context.stroke();
    context.closePath();
    
}
function draw_wall(canvas_id){
// 获得 canvas.context
clearCanvas(canvas_id)
get_design_data()
polygon(canvas_id,submit_data["design"]["coords"],);
}
function add_tu(){
    //展示数据
    let html = `<tr>
        <td>${document.querySelector("#po_h").value}</td>
        <td>${document.querySelector("#po_v").value}</td>
        <td>${document.querySelector("#tu_h").value}</td>
        <td>${document.querySelector("#tu_w").value}</td>
        <td>${document.querySelector("#tu_j").value}</td>
        <td><input type="button" class="btn btn-danger" value="删除此行" onclick="delate(this); draw_po()"></td>
    </tr>`
    document.querySelector("#show_po").innerHTML += html
    document.querySelector("#show_po_div").style.display = "block"
}
function delate(e){
e.parentElement.parentElement.innerHTML=""
}
function get_y(x0,y0,x1,y1,x){
    //计算y值
    let k = (y1-y0)/(x1-x0)
    let y = k*(x-x0)+y0
    return y
}
function get_tu_data(){
    let tu = [submit_data["design"]["coord_qd"]]
    let tuzhu = []
    let trs = document.querySelectorAll("#show_po > tr")
    for(let i=0;i<trs.length;i++){
        let tds = trs[i].querySelectorAll("td")
        if((! tds) ||(!tds[0]) ){
            continue
        }
        if (tds[0].innerText && tds[1].innerText ){
            let x = parseFloat(tds[0].innerText)+tu[tu.length-1][0]
            let y = parseFloat(tds[1].innerText)+tu[tu.length-1][1]
            tu.push([x,y])
        }
        if(tds[0].innerText && tds[1].innerText &&tds[2].innerText && tds[3].innerText && tds[4].innerText){
            let height = parseFloat(tds[2].innerText)
            let width = parseFloat(tds[3].innerText)
            let j = parseFloat(tds[4].innerText)
            let x = j + tu[tu.length-2][0]
            let y = get_y(tu[tu.length-2][0],tu[tu.length-2][1],tu[tu.length-1][0],tu[tu.length-1][1],x)
            let x1 = x + width
            let y1 = get_y(tu[tu.length-2][0],tu[tu.length-2][1],tu[tu.length-1][0],tu[tu.length-1][1],x1)
            tuzhu.push({
                "start":[x,y],
                "end":[x1,y1],
                "height":height,
                "width":width,
                "coords":[[x,y],[x,y+height],[x1,y1+height],[x1,y1],[x,y]]
            })
        }
    }
    submit_data["tu"] = {
        "po":tu,
        "tuzhu":tuzhu
    }
    return submit_data["tu"]
}
function draw_po(canvas_id){
    draw_wall(canvas_id)
    let tu = get_tu_data()
   
    if(tu["po"].length>1){
        polygon(canvas_id,tu["po"])
    }
    if(tu["tuzhu"].length>0){
        for(let i=0;i<tu["tuzhu"].length;i++){
            let tuzhu = tu["tuzhu"][i]
            let x = tuzhu["start"][0]
            let y = tuzhu["start"][1]
            let x1 = tuzhu["end"][0]
            let y1 = tuzhu["end"][1]
            let height = tuzhu["height"]
            polygon(canvas_id,[[x,y],[x,y+height],[x1,y1+height],[x1,y1],[x,y]])
        }
    }
    // console.log(submit_data)
}

function get_parms(){
    let parms = {}
    for(i of ["accuracy","f","gamma_wall","base_type","base_sigma","phi","gamma_soil","C_z","k_H","gamma_0"]){
        if(document.getElementById(i)&& document.getElementById(i).value){
            if(i=="base_type"){
                parms[i] = document.getElementById(i).value
            }else{
                parms[i] = parseFloat(document.getElementById(i).value)
            }
        }else{
            console.log("未输入",i)
            alert("请输入 "+document.getElementById(i).parentElement.innerText)
            return false
        }
    }
    submit_data["parms"] = parms
    return parms
}

function write_check(check,data){
    for(i of ["K_c","K_0","e","sigma_1k","sigma_2k","sigma_pk"]){
        let load =/sigma/.test(i)? data[check]["sigma"][i] : data[check][i]
        let ok_str = `<span><i>${safe}</i> 通过</span>`
        if(load[1] < load[0]){
            ok_str = `<span> <i>${not_safe}</i>不通过</span>`
            check_flag  = false
        }
        $("#"+check).find("#result_"+i).append( `
     <td  class="result_table" >${load[0]}</td>
     <td  class="result_table" >${load[1]}</td>
     <td  class="result_table" >${ok_str}</td>`)
    }
}
function submit_data_to_server(){
    //用ajax将submit上传到服务器
    get_design_data()
    get_tu_data()
    let parms = get_parms()
    if(!parms){
        return
    }
    if(submit_data["tu"]["po"].length<2){
        alert("请输入坡面数据")
        return
    }
    $("#mask").show()
    // console.log("submit_data",JSON.stringify(submit_data))
    $.ajax({
        url: "/submit",
        type: "POST",
        data: JSON.stringify(submit_data),
        contentType: "application/json; charset=utf-8",
        dataType: "json",
        success: function(data){
            $("#mask").hide()
            window.data = data
            if(data["code"]!=0){
                alert(data["msg"])
                return
            }
            $(".result_text").empty()
            $(".result_table").remove()
            $("#result_show").show()
            data = data["data"]
            polygon("quad",[submit_data["design"]["coord_qz"],data["jx_start"]]);
            polygon("quad",[submit_data["design"]["coord_qz"],data["jx_end"]]);
            polygon("quad",[submit_data["design"]["coord_qz"],submit_data["design"]["coord_qd"]]);
            $("#stress_chart").attr("href",'.'+data["png_stress"]);

            fig =  JSON.parse(data["png_E_3d"])
            Plotly.newPlot("result3d", fig)

            for(i of ["theta","alpha","B","area_G","area_W","area_W1","sum_M_y","sum_M_0","sum_M_0","sum_N"]){
                $("#result_"+i).html(`${data[i]}`)
            }
            if(!data["have_2_plm"]) {
                $("#result_have_2_plm").html(`假想墙背夹角：`)
                $("#result_area_W_desc").html(`真实墙背与假想墙背土体截面积：`)
                $("#result_area_W1_desc").html(`第一破裂面与假想墙背土体截面积：`)
            }else{
                $("#result_have_2_plm").html(`第二破裂面夹角：`)
                $("#result_area_W_desc").html(`真实墙背与第二破裂面土体截面积：`)
                $("#result_area_W1_desc").html(`第一破裂面与第二破裂面土体截面积：`)
            }
            for (i of ["E_a","G","W_toe","W","W1","G_p","W_p","W1_p"]){
            $("#result_"+i).append( `<td class="result_table" >${data[i].num}</td>
                                     <td class="result_table" >${data[i].x}</td>
                                     <td class="result_table" >${data[i].y}</td>
                                     <td class="result_table" >${(String(data[i].pos))}</td>`)
            }
            check_flag = true
            for(i of ["check_1","check_4"]) write_check(i,data)
            if (check_flag){
                $("#result_con").html("<span class='success'>经过计算，该挡墙符合相关设计规范要求，稳定性验算通过</span>")
            }else{
                $("#result_con").html("<span class = 'danger'>经过计算，该挡墙不符合相关设计规范要求，稳定性验算不通过</span>")
            }
            MathJax.Hub.Queue(["Typeset",MathJax.Hub]);
            $("html,body").animate({scrollTop: $("#result_show").offset().top},300);

        },
        error: function(data){
            $("#mask").hide()
            console.log(data)
            alert("程序错误！！")
        }

    });
}