import os

_directory_stack =[]

def pushd(target,mkdir=False):
    oldpath = os.getcwd()
    newpath = os.path.abspath(target) 

    _directory_stack.append(oldpath)
    if not os.path.exists(newpath) and mkdir==True:
        print("Making directory: %s" % (newpath))
        os.mkdir(newpath)

    if os.path.exists(newpath) :
        os.chdir(newpath)
    else:
        print("Cant cd to new path %s" % (target))
    print("Current directory: %s" % (os.getcwd()) )

def popd():
    try: 
        current_top = _directory_stack.pop()
        print("Changing directory to %s" % (os.getcwd()) )
        os.chdir(current_top)
    except:
        print("Directory stack empty or directory does not exist")

if __name__ == "__main__":
    print("hworld") 
    pushd("./aatest") 
    pushd("./aatest", mkdir=True) 
    pushd("./aatest") 
    popd() 
    popd() 
    popd()




