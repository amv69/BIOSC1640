from flask import Flask, render_template
app = Flask(__name__)

data = [
    {
        'fastA': 'ATCG',
        'lengthOne': '5',
        'lengthTwo': '10',
        'tempOne': '20',
        'tempTwo': '40',
    }

]

@app.route("/")
@app.route("/home")
def hello():
    return render_template('home.html', data=data)

@app.route("/about")
def about():
    return render_template('about.html', title='About')

if __name__ == '__main__':
    app.run(debug=True)


