from starlette.responses import RedirectResponse
from fastapi import Body, FastAPI
from fastapi.staticfiles import StaticFiles
from handle import *

app = FastAPI()
app.mount("/static", StaticFiles(directory="./html"), name="static")  
@app.get("/index.html")
async def redirect_typer():
    return RedirectResponse("/static/index.html")
@app.get("/")
async def redirect_typer():
    return RedirectResponse("/static/index.html")

@app.post('/submit')
def index(data=Body(...)):
	try:
		tmp = Dangtuqiang(data)
		return_data = tmp.get_result()
		return {
			"code":0,
			"data":return_data,
			"msg":"success"
		}
	except Exception as e:
		return{
			"code":1,
			"msg":"程序错误！\n"+str(e)
		}
	
    
 